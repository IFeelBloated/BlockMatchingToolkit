import vapoursynth as vs
import mvmulti
import math

### Global Settings ###
fmtc_args         = dict (fulls=True, fulld=True)
msuper_args       = dict (chroma=False, hpad=32, vpad=32, sharp=2, levels=0)
manalyze_args     = dict (search=3, chroma=False, truemotion=True, trymany=True, levels=0, badrange=-24)
mrecalculate_args = dict (chroma=False, truemotion=True, search=3, smooth=1)
mdegrain_args     = dict (plane=0, limit=1.0)
nnedi_args        = dict (field=1, dh=True, nns=4, qual=2, etype=1, nsize=0)
canny_args        = dict (sigma=1.5, mode=1, op=0)
deconv_args       = dict (line=0, wn=0.99, x=1, y=1, fr=25, scale=0.0059)

### Helpers ###
def gauss (src, p=30):
    core          = vs.get_core ()
    resample      = core.fmtc.resample
    upsmp         = resample (src, src.width * 2, src.height * 2, kernel="gauss", a1=100, **fmtc_args)
    clip          = resample (upsmp, src.width, src.height, kernel="gauss", a1=p, **fmtc_args)
    return clip

def padding (src, left=0, right=0, top=0, bottom=0):
    core          = vs.get_core ()
    resample      = core.fmtc.resample
    w             = src.width
    h             = src.height
    clip          = resample (src, w+left+right, h+top+bottom, -left, -top, w+left+right, h+top+bottom, kernel="point", **fmtc_args)
    return clip

def xymax (src1, src2):
    core          = vs.get_core ()
    Expr          = core.std.Expr
    clip          = Expr ([src1, src2], ["x y max"])
    return clip

def xymin (src1, src2):
    core          = vs.get_core ()
    Expr          = core.std.Expr
    clip          = Expr ([src1, src2], ["x y min"])
    return clip

def max_dif (src1, src2, ref):
    core          = vs.get_core ()
    Expr          = core.std.Expr
    clip          = Expr ([src1, src2, ref], ["x z - abs y z - abs > x y ?"])
    return clip

def min_dif (src1, src2, ref):
    core          = vs.get_core ()
    Expr          = core.std.Expr
    clip          = Expr ([src1, src2, ref], ["x z - abs y z - abs > y x ?"])
    return clip

def thr_merge (flt, src, ref=None, thr=0.0009765625, elast=None):
    core          = vs.get_core ()
    Expr          = core.std.Expr
    MakeDiff      = core.std.MakeDiff
    MergeDiff     = core.std.MergeDiff
    ref           = src if ref is None else ref
    elast         = thr / 2 if elast is None else elast
    BExp          = ["x {thr} {elast} + z - 2 {elast} * / * y {elast} z + {thr} - 2 {elast} * / * +".format (thr=thr, elast=elast)]
    BDif          = Expr (src, "0.0")
    PDif          = Expr ([flt, src], "x y - 0.0 max")
    PRef          = Expr ([flt, ref], "x y - 0.0 max")
    PBLD          = Expr ([PDif, BDif, PRef], BExp)
    NDif          = Expr ([flt, src], "y x - 0.0 max")
    NRef          = Expr ([flt, ref], "y x - 0.0 max")
    NBLD          = Expr ([NDif, BDif, NRef], BExp)
    BLDD          = MakeDiff (PBLD, NBLD)
    BLD           = MergeDiff (src, BLDD)
    UDN           = Expr ([flt, ref, BLD], ["x y - abs {thr} {elast} - > z x ?".format (thr=thr, elast=elast)])
    clip          = Expr ([flt, ref, UDN, src], ["x y - abs {thr} {elast} + < z a ?".format (thr=thr, elast=elast)])
    return clip

def clamp (src, bright_limit, dark_limit, overshoot=0.0, undershoot=0.0):
    core          = vs.get_core ()
    Expr          = core.std.Expr
    clip          = Expr ([src, bright_limit, dark_limit], ["x y {os} + > y {os} + x ? z {us} - < z {us} - x ?".format (os=overshoot, us=undershoot)])
    return clip

def maxmulti (src, start=None, a=0, tr=6):
    core          = vs.get_core ()
    Expr          = core.std.Expr
    SelectEvery   = core.std.SelectEvery
    start         = Expr (SelectEvery (src, tr*2+1, 0), ["0.0"]) if start is None else start
    max           = xymax (start, SelectEvery (src, tr*2+1, a))
    a             = a+1
    clip          = max if a == tr*2+1 else maxmulti (src, start=max, tr=tr, a=a)
    return clip

def minmulti (src, start=None, a=0, tr=6):
    core          = vs.get_core ()
    Expr          = core.std.Expr
    SelectEvery   = core.std.SelectEvery
    start         = Expr (SelectEvery (src, tr*2+1, 0), ["1.0"]) if start is None else start
    min           = xymin (start, SelectEvery (src, tr*2+1, a))
    a             = a+1
    clip          = min if a == tr*2+1 else minmulti (src, start=min, tr=tr, a=a)
    return clip

def hipass (src, sharp, p=16):
    core          = vs.get_core ()
    MakeDiff      = core.std.MakeDiff
    MergeDiff     = core.std.MergeDiff
    hif           = MakeDiff (sharp, gauss (sharp, p=p))
    clip          = MergeDiff (gauss (src, p=p), hif)
    return clip

### Denoisers ###
def halonr (src, rclip=None, a=32, h=6.4, thr=0.00390625, elast=None):
    core          = vs.get_core ()
    NLMeans       = core.knlm.KNLMeansCL
    Canny         = core.tcanny.TCanny
    Expand        = core.flt.Maximum
    Inflate       = core.flt.Inflate
    Crop          = core.std.CropRel
    Expr          = core.std.Expr
    MaskedMerge   = core.std.MaskedMerge
    elast         = thr / 6 if elast is None else elast
    pad           = padding (src, a, a, a, a)
    padref        = pad if rclip is None else padding (rclip, a, a, a, a)
    Clean         = NLMeans (pad, d=0, a=a, s=0, h=h, rclip=padref)
    CLNRef        = 0 if rclip is None else NLMeans (padref, d=0, a=a, s=0, h=h)
    EMask         = Canny (Clean, **canny_args) if rclip is None else Canny (CLNRef, **canny_args)
    EMask         = Expr (EMask, "x 0.24 - 3.2 * 0.0 max 1.0 min")
    EMask         = Expand (EMask)
    EMask         = Inflate (EMask)
    MRG           = MaskedMerge (pad, Clean, EMask)
    Final         = thr_merge (pad, MRG, thr=thr, elast=elast)
    clip          = Crop (Final, a, a, a, a)
    return clip

def reconstructor (src, pelclip=None, nrlevel=1, a=32, h1=None, h2=None, divide=4, tr=6, pel=4, dct=5, thsad=2000, thscd1=10000, thscd2=255):
    core          = vs.get_core ()
    NLMeans       = core.knlm.KNLMeansCL
    Median        = core.flt.Median
    Repair        = core.rgsf.Repair
    MakeDiff      = core.std.MakeDiff
    MergeDiff     = core.std.MergeDiff
    Crop          = core.std.CropRel
    Expr          = core.std.Expr
    MSuper        = core.mvsf.Super
    MAnalyze      = mvmulti.Analyze
    MDegrainN     = mvmulti.DegrainN
    if (nrlevel == 1):
       h1         = 9.6 if h1 is None else h1
       h2         = 6.4 if h2 is None else h2
    if (nrlevel == 2):
       h1         = 128.0 if h1 is None else h1
       h2         = 7.2 if h2 is None else h2
    hfine         = h2 * 2 / divide
    den           = divide - 1
    windmin       = a // pow (2, den)
    supsrh        = MSuper (src, pelclip=pelclip, rfilter=4, pel=pel, **msuper_args)
    suprdr        = MSuper (src, pelclip=pelclip, rfilter=2, pel=pel, **msuper_args)
    vmulti        = MAnalyze (supsrh, overlap=2, blksize=4, divide=0, tr=tr, dct=dct, **manalyze_args)
    clean         = MDegrainN (src, suprdr, vmulti, tr=tr, thsad=thsad, thscd1=thscd1, thscd2=thscd2, **mdegrain_args)
    ref           = padding (Repair (clean, src, mode=17), a+1, a+1, a+1, a+1)
    dif           = padding (MakeDiff (src,clean), a+1, a+1, a+1, a+1)
    difcln        = NLMeans (dif, d=0, a=windmin, s=1, h=h1, rclip=ref)
    difcln        = gauss (Median (difcln), p=16)
    dif           = Crop (Expr ([dif, difcln], "x 0.5 - y 0.5 - * 0 < 0.5 x 0.5 - abs y 0.5 - abs < x y ? ?"), a+1, a+1, a+1, a+1)
    mrg           = MergeDiff (clean, dif)
    rep           = Repair (mrg, src, mode=17)
    dif           = padding (MakeDiff (src, rep), a+1, a+1, a+1, a+1)
    difcln        = NLMeans (dif, d=0, a=windmin, s=1, h=h1, rclip=padding (rep, a+1, a+1, a+1, a+1))
    difcln        = Crop (difcln, a+1, a+1, a+1, a+1)
    cln           = MergeDiff (rep, difcln)
    pad           = padding (src, a+1, a+1, a+1, a+1) if nrlevel == 1 else padding (cln, a+1, a+1, a+1, a+1)
    init          = pad if nrlevel == 2 else padding (cln, a+1, a+1, a+1, a+1)
    init          = NLMeans (init, d=0, a=windmin, s=1, h=h2) if nrlevel == 2 else NLMeans (init, d=0, a=windmin, s=1, h=hfine)
    def inline (flt, n):
        str       = n * h2 / den + hfine * (1 - n / den)
        window    = a // pow (2, den - n)
        flt       = init if n == den else flt
        dif       = MakeDiff (pad, flt)
        dif       = NLMeans (dif, d=0, a=window, s=1, h=str, rclip=flt)
        fnl       = MergeDiff (flt, dif)
        n         = n - 1
        return fnl if n == -1 else inline (fnl, n)
    Final         = inline (None, den)
    clip          = Crop (Final, a+1, a+1, a+1, a+1)
    return clip

def generalnr (src, srclow=None, a=32, h=1.2, sigma=8.0, block_size=8, block_step=1, group_size=32, bm_range=24, bm_step=1, lowpass=8):
    core          = vs.get_core ()
    NLMeans       = core.knlm.KNLMeansCL
    BM3DBasic     = core.bm3d.Basic
    BM3DFinal     = core.bm3d.Final
    Crop          = core.std.CropRel
    srclow        = src if srclow is None else srclow
    pad           = padding (src, a+4, a+4, a+4, a+4)
    hflt          = NLMeans (pad, d=0, a=a, s=4, h=h)
    hif           = Crop (hflt, a+4, a+4, a+4, a+4)
    ref           = BM3DBasic (srclow, sigma=sigma, block_size=block_size, block_step=block_step, group_size=group_size, bm_range=bm_range, bm_step=bm_step)
    lowf          = BM3DFinal (srclow, ref, sigma=sigma, block_size=block_size, block_step=block_step, group_size=group_size, bm_range=bm_range, bm_step=bm_step)
    clip          = hipass (lowf, hif, p=lowpass)
    return clip

def nrfinal (spatial, dif, vmulti, peldif=None, pel=4, tr=6, thsad=10000, thscd1=10000, thscd2=255, repmode=13):
    core          = vs.get_core ()
    Repair        = core.rgsf.Repair
    Expr          = core.std.Expr
    MergeDiff     = core.std.MergeDiff
    MSuper        = core.mvsf.Super
    MDegrainN     = mvmulti.DegrainN
    superclip     = MSuper (dif, pelclip=peldif, rfilter=2, pel=pel, **msuper_args)
    blankd        = Expr ([dif], "0.5")
    mc            = MDegrainN (blankd, superclip, vmulti, tr=tr, thsad=thsad, thscd1=thscd1, thscd2=thscd2, **mdegrain_args)
    nr            = MergeDiff (spatial, mc)
    repclip       = MergeDiff (spatial, dif)
    clip          = Repair (nr, repclip, repmode)
    return clip

### Sharpeners ###
def delicatesharp (src, relaxed=False):
    core          = vs.get_core ()
    Median        = core.flt.Median
    Repair        = core.rgsf.Repair
    RemoveGrain   = core.rgsf.RemoveGrain
    Expr          = core.std.Expr
    MakeDiff      = core.std.MakeDiff
    MergeDiff     = core.std.MergeDiff
    def inline (src):
        rg11      = RemoveGrain (src, 11)
        rg11D     = MakeDiff (src, rg11)
        rg11DR    = RemoveGrain (rg11D, 11)
        rg11DD    = Expr ([rg11D, rg11DR], ["x y - x 0.5 - * 0 < 0.5 x y - abs x 0.5 - abs < x y - 0.5 + x ? ?"])
        clip      = MakeDiff (src, rg11DD)
        return clip
    blur          = Repair (Median (src), src, mode=16) if relaxed else Median (src)
    dif           = Expr ([src, blur], ["y x - 3 * 0.5 +"])
    dif           = inline (dif)
    difrg         = RemoveGrain (dif, 11)
    fdif          = Expr ([dif, difrg], ["y 0.5 - abs x 0.5 - abs > y 0.5 ?"])
    clip          = MergeDiff (src, fdif)
    return clip

def nrsharp (src):
    core          = vs.get_core ()
    NNEDI         = core.nnedi3.nnedi3
    Median        = core.flt.Median
    resample      = core.fmtc.resample
    MakeDiff      = core.std.MakeDiff
    MergeDiff     = core.std.MergeDiff
    Transpose     = core.std.Transpose
    w             = src.width
    h             = src.height
    def inline (clp, n):
        blur      = resample (clp, w*8, h*8, kernel="cubic", a1=1, a2=0, **fmtc_args)
        sharp     = resample (clp, w*8, h*8, kernel="cubic", a1=-1, a2=0, **fmtc_args)
        dif       = MakeDiff (blur, sharp)
        dif       = resample (dif, w*2, h*2, kernel="cubic", a1=-1, a2=0, **fmtc_args)
        clp       = MergeDiff (clp, dif)
        n         = n-1
        return clp if n == 0 else inline (clp, n)
    u2x           = Transpose (NNEDI (Transpose (NNEDI (src, **nnedi_args)), **nnedi_args))
    ublur         = inline (u2x, 4)
    cdif          = resample (MakeDiff (u2x, ublur), w, h, -0.25, -0.25, kernel="cubic", a1=-4, a2=-2, **fmtc_args)
    csharp        = hipass (src, MergeDiff (src, cdif), p=16)
    msharp        = MergeDiff (src, MakeDiff (src, Median (src)))
    clip          = min_dif (csharp, msharp, src)
    return clip

def deconvolution (src):
    core          = vs.get_core ()
    Deconv        = core.vcfreq.Sharp
    dcv           = Deconv (Deconv (src, **deconv_args), **deconv_args)
    clip          = hipass (src, dcv, p=8)
    return clip

def sharpfinal (soft, dif, limit, vmulti, peldif=None, pellimit=None, pel=4, tr=6, thsadA=10000, thsadL=400, thscd1=10000, thscd2=255, a=32, h=6.4, thr=0.00390625, elast=None, repmode=13, str=1.00):
    core          = vs.get_core ()
    Repair        = core.rgsf.Repair
    Maximum       = core.flt.Maximum
    Minimum       = core.flt.Minimum
    Expr          = core.std.Expr
    MakeDiff      = core.std.MakeDiff
    MergeDiff     = core.std.MergeDiff
    MSuper        = core.mvsf.Super
    MDegrainN     = mvmulti.DegrainN
    MCompensate   = mvmulti.Compensate
    expression    = ["{x} {y} {x} - abs 4 / 1 4 / pow 4 * 1.6 * {str} * {y} {x} - {y} {x} - abs 1.001 + / * + 256 /".format (str=str, x="x 256 *", y="y 256 *")]
    blankd        = Expr ([dif], ["0.5"])
    superdif      = MSuper (dif, pelclip=peldif, rfilter=2, pel=pel, **msuper_args)
    supercmp      = MSuper (limit, pelclip=pellimit, rfilter=2, pel=pel, **msuper_args)
    MDG           = MDegrainN (blankd, superdif, vmulti, tr=tr, thsad=thsadA, thscd1=thscd1, thscd2=thscd2, **mdegrain_args)
    TA            = Repair (MergeDiff (soft, MDG), MergeDiff (soft, dif), mode=repmode)
    comp          = MCompensate (soft, supercmp, vmulti, tr=tr, thsad=thsadL, thscd1=thscd1, thscd2=thscd2)
    max           = maxmulti (comp, tr=tr)
    min           = minmulti (comp, tr=tr)
    TL            = clamp (TA, max, min, overshoot=0.0, undershoot=0.0)
    local         = Expr ([soft, TL], expression)
    bright        = Maximum (limit)
    dark          = Minimum (limit)
    SL            = clamp (local, bright, dark, overshoot=0.0, undershoot=0.0)
    SRPD          = halonr (MakeDiff (SL, soft), rclip=SL, a=a, h=h, thr=thr, elast=elast)
    clip          = MergeDiff (soft, SRPD)
    return clip

### Motion Estimation ###
def getvectors (src, pelclip=None, tr=6, pel=4, dct=5, thsad=400):
    core          = vs.get_core ()
    MSuper        = core.mvsf.Super
    MAnalyze      = mvmulti.Analyze
    MRecalculate  = mvmulti.Recalculate
    supersoft     = MSuper (src, pelclip=pelclip, rfilter=4, pel=pel, **msuper_args)
    supersharp    = MSuper (src, pelclip=pelclip, rfilter=2, pel=pel, **msuper_args)
    vmulti        = MAnalyze (supersoft, overlap=16, blksize=32, divide=2, tr=tr, dct=dct, **manalyze_args)
    vmulti        = MRecalculate (supersoft, vmulti, overlap=8, blksize=16, divide=2, tr=tr, thsad=thsad/2, dct=dct, **mrecalculate_args)
    vmulti        = MRecalculate (supersharp, vmulti, overlap=4, blksize=8, divide=2, tr=tr, thsad=thsad/2, dct=dct, **mrecalculate_args)
    vmulti        = MRecalculate (supersharp, vmulti, overlap=2, blksize=4, divide=0, tr=tr, thsad=thsad/2, dct=dct, **mrecalculate_args)
    return vmulti

def genpelclip (src, pel=4):
    core          = vs.get_core ()
    NNEDI         = core.nnedi3.nnedi3
    Transpose     = core.std.Transpose
    u2x           = Transpose (NNEDI (Transpose (NNEDI (src, **nnedi_args)), **nnedi_args))
    u4x           = Transpose (NNEDI (Transpose (NNEDI (u2x, **nnedi_args)), **nnedi_args))
    clip          = u2x if pel == 2 else u4x
    return clip
