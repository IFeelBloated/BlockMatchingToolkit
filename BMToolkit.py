import vapoursynth as vs
import mvmulti
import math

### Global Settings ###
fmtc_args           = dict (fulls=True, fulld=True)
msuper_args         = dict (chroma=False, hpad=32, vpad=32, sharp=2, levels=0)
manalyze_args       = dict (search=3, chroma=False, truemotion=True, trymany=True, levels=0, badrange=-24)
mrecalculate_args   = dict (chroma=False, truemotion=True, search=3, smooth=1)
mdegrain_args       = dict (plane=0, limit=1.0)
nnedi_args          = dict (field=1, dh=True, nns=4, qual=2, etype=1, nsize=0)
canny_args          = dict (sigma=1.5, mode=1, op=0)
deconv_args         = dict (line=0, wn=0.99, x=1, y=1, fr=25, scale=0.0059)

### Helpers ###
def gauss (src, p=30):
    core            = vs.get_core ()
    resample        = core.fmtc.resample
    upsmp           = resample (src, src.width * 2, src.height * 2, kernel="gauss", a1=100, **fmtc_args)
    clip            = resample (upsmp, src.width, src.height, kernel="gauss", a1=p, **fmtc_args)
    return clip

def padding (src, left=0, right=0, top=0, bottom=0):
    core            = vs.get_core ()
    resample        = core.fmtc.resample
    w               = src.width
    h               = src.height
    clip            = resample (src, w+left+right, h+top+bottom, -left, -top, w+left+right, h+top+bottom, kernel="point", **fmtc_args)
    return clip

def xymax (src1, src2):
    core            = vs.get_core ()
    Expr            = core.std.Expr
    clip            = Expr ([src1, src2], ["x y max"])
    return clip

def xymin (src1, src2):
    core            = vs.get_core ()
    Expr            = core.std.Expr
    clip            = Expr ([src1, src2], ["x y min"])
    return clip

def max_dif (src1, src2, ref):
    core            = vs.get_core ()
    Expr            = core.std.Expr
    clip            = Expr ([src1, src2, ref], ["x z - abs y z - abs > x y ?"])
    return clip

def min_dif (src1, src2, ref):
    core            = vs.get_core ()
    Expr            = core.std.Expr
    clip            = Expr ([src1, src2, ref], ["x z - abs y z - abs > y x ?"])
    return clip

def thr_merge (flt, src, ref=None, thr=0.0009765625, elast=None):
    core            = vs.get_core ()
    Expr            = core.std.Expr
    MakeDiff        = core.std.MakeDiff
    MergeDiff       = core.std.MergeDiff
    ref             = src if ref is None else ref
    elast           = thr / 2 if elast is None else elast
    BExp            = ["x {thr} {elast} + z - 2 {elast} * / * y {elast} z + {thr} - 2 {elast} * / * +".format (thr=thr, elast=elast)]
    BDif            = Expr (src, "0.0")
    PDif            = Expr ([flt, src], "x y - 0.0 max")
    PRef            = Expr ([flt, ref], "x y - 0.0 max")
    PBLD            = Expr ([PDif, BDif, PRef], BExp)
    NDif            = Expr ([flt, src], "y x - 0.0 max")
    NRef            = Expr ([flt, ref], "y x - 0.0 max")
    NBLD            = Expr ([NDif, BDif, NRef], BExp)
    BLDD            = MakeDiff (PBLD, NBLD)
    BLD             = MergeDiff (src, BLDD)
    UDN             = Expr ([flt, ref, BLD], ["x y - abs {thr} {elast} - > z x ?".format (thr=thr, elast=elast)])
    clip            = Expr ([flt, ref, UDN, src], ["x y - abs {thr} {elast} + < z a ?".format (thr=thr, elast=elast)])
    return clip

def clamp (src, bright_limit, dark_limit, overshoot=0.0, undershoot=0.0):
    core            = vs.get_core ()
    Expr            = core.std.Expr
    clip            = Expr ([src, bright_limit, dark_limit], ["x y {os} + > y {os} + x ? z {us} - < z {us} - x ?".format (os=overshoot, us=undershoot)])
    return clip

def maxmulti (src, start=None, a=0, tr=6):
    core            = vs.get_core ()
    Expr            = core.std.Expr
    SelectEvery     = core.std.SelectEvery
    start           = Expr (SelectEvery (src, tr*2+1, 0), ["0.0"]) if start is None else start
    max             = xymax (start, SelectEvery (src, tr*2+1, a))
    a               = a+1
    clip            = max if a == tr*2+1 else maxmulti (src, start=max, tr=tr, a=a)
    return clip

def minmulti (src, start=None, a=0, tr=6):
    core            = vs.get_core ()
    Expr            = core.std.Expr
    SelectEvery     = core.std.SelectEvery
    start           = Expr (SelectEvery (src, tr*2+1, 0), ["1.0"]) if start is None else start
    min             = xymin (start, SelectEvery (src, tr*2+1, a))
    a               = a+1
    clip            = min if a == tr*2+1 else minmulti (src, start=min, tr=tr, a=a)
    return clip

def hipass (src, sharp, p=16):
    core            = vs.get_core ()
    MakeDiff        = core.std.MakeDiff
    MergeDiff       = core.std.MergeDiff
    hif             = MakeDiff (sharp, gauss (sharp, p=p))
    clip            = MergeDiff (gauss (src, p=p), hif)
    return clip

def genblockmask (src):
    core            = vs.get_core ()
    resample        = core.fmtc.resample
    Expr            = core.std.Expr
    BlankClip       = core.std.BlankClip
    AddBorders      = core.std.AddBorders
    StackHorizontal = core.std.StackHorizontal
    StackVertical   = core.std.StackVertical
    CropAbs         = core.std.CropAbs
    clip            = BlankClip (src, 24, 24, color=0.0)
    clip            = AddBorders (clip, 4, 4, 4, 4, color=1.0)
    clip            = StackHorizontal ([clip, clip, clip, clip])
    clip            = StackVertical ([clip, clip, clip, clip])
    clip            = resample (clip, 32, 32, kernel="point", **fmtc_args)
    clip            = Expr (clip, ["x 0.0 > 1.0 0.0 ?"])
    clip            = StackHorizontal ([clip, clip, clip, clip, clip, clip, clip, clip])
    clip            = StackVertical ([clip, clip, clip, clip, clip, clip])
    clip            = StackHorizontal ([clip, clip, clip, clip, clip, clip])
    clip            = StackVertical ([clip, clip, clip, clip, clip])
    clip            = CropAbs (clip, src.width, src.height, 0, 0)
    return clip

### Denoisers ###
def halonr (src, a=32, h=6.4, thr=0.00390625, elast=None, lowpass=8):
    core            = vs.get_core ()
    NLMeans         = core.knlm.KNLMeansCL
    Canny           = core.tcanny.TCanny
    Expand          = core.flt.Maximum
    Inflate         = core.flt.Inflate
    Crop            = core.std.CropRel
    Expr            = core.std.Expr
    MaskedMerge     = core.std.MaskedMerge
    elast           = thr / 6 if elast is None else elast
    pad             = padding (src, a, a, a, a)
    Clean           = hipass (pad, NLMeans (pad, d=0, a=a, s=0, h=h, wref=1.0, rclip=pad), p=lowpass)
    EMask           = Canny (Clean, **canny_args)
    EMask           = Expr (EMask, "x 0.24 - 128.0 * 0.0 max 1.0 min")
    EMask           = Expand (Inflate (Expand (EMask)))
    MRG             = MaskedMerge (pad, Clean, EMask)
    Final           = thr_merge (pad, MRG, thr=thr, elast=elast)
    clip            = Crop (Final, a, a, a, a)
    return clip

def crapnr (src, pelclip=None, nrlevel=1, deblock=True, h=6.4, thr=0.03125, elast=0.015625, tr=6, pel=4, thsad=2000, thscd1=10000, thscd2=255, sigma=8.0, block_size=8, block_step=1, group_size=32, bm_range=24, bm_step=1):
    core            = vs.get_core ()
    NLMeans         = core.knlm.KNLMeansCL
    BM3DBasic       = core.bm3d.Basic
    BM3DFinal       = core.bm3d.Final
    MakeDiff        = core.std.MakeDiff
    MergeDiff       = core.std.MergeDiff
    MaskedMerge     = core.std.MaskedMerge
    Crop            = core.std.CropRel
    MSuper          = core.mvsf.Super
    MAnalyze        = mvmulti.Analyze
    MDegrainN       = mvmulti.DegrainN
    hfine           = pow (1.1988568728336214663622280225868, h)
    def inline_BM_inter ():
        supsrh      = MSuper (src, pelclip=pelclip, rfilter=4, pel=pel, **msuper_args)
        suprdr      = MSuper (src, pelclip=pelclip, rfilter=2, pel=pel, **msuper_args)
        vmulti      = MAnalyze (supsrh, overlap=2, blksize=4, divide=0, tr=tr, dct=5, **manalyze_args)
        return MDegrainN (src, suprdr, vmulti, tr=tr, thsad=thsad, thscd1=thscd1, thscd2=thscd2, **mdegrain_args)
    def inline_BM_intra (src, h):
        ref         = BM3DBasic (src, sigma=sigma, block_size=block_size, block_step=block_step, group_size=group_size, bm_range=bm_range, bm_step=bm_step) 
        bm3d        = BM3DFinal (src, ref, sigma=sigma, block_size=block_size, block_step=block_step, group_size=group_size, bm_range=bm_range, bm_step=bm_step)
        return NLMeans (bm3d, d=0, a=block_size//2, s=1, h=h, wref=1.0)
    def inline_NLM (flt, init, src, n):
        c1          = 1.0707892518365290738330599429051
        c2          = 0.4798695862246764421520306169363
        str         = n * h / 4 + hfine * (1 - n / 4)
        weight      = pow (c1, str * (4 - n)) - c2
        window      = 32 // pow (2, n)
        flt         = init if n == 4 else flt
        dif         = MakeDiff (src, flt)
        dif         = NLMeans (dif, d=0, a=window, s=1, h=str, wref=weight, rclip=flt)
        fnl         = MergeDiff (flt, dif)
        n           = n - 1
        return fnl if n == -1 else inline_NLM (fnl, init, src, n)
    pad             = padding (src, 33, 33, 33, 33)
    BM_inter        = inline_BM_inter ()
    BM_inter_fine   = Crop (inline_NLM (None, padding (BM_inter, 33, 33, 33, 33), pad, 4), 33, 33, 33, 33) if nrlevel == 1 else 0
    BM_inter_raw    = thr_merge (BM_inter_fine, hipass (BM_inter_fine, BM_inter, 8), thr=thr, elast=elast) if nrlevel == 1 and deblock else 0
    BM_inter        = MaskedMerge (BM_inter_fine, BM_inter_raw, genblockmask (src)) if nrlevel == 1 and deblock else BM_inter_fine if nrlevel == 1 else BM_inter
    BM_inter        = padding (BM_inter, 33, 33, 33, 33)
    BM_intra        = inline_BM_intra (BM_inter, hfine) if nrlevel == 1 else inline_BM_intra (BM_inter, h)
    BM_intra        = inline_NLM (None, BM_intra, BM_inter, 4)
    clip            = Crop (BM_intra, 33, 33, 33, 33)
    return clip

def generalnr (src, srclow=None, a=32, h=2.4, sigma=8.0, block_size=8, block_step=1, group_size=32, bm_range=24, bm_step=1, lowpass=8):
    core            = vs.get_core ()
    NLMeans         = core.knlm.KNLMeansCL
    BM3DBasic       = core.bm3d.Basic
    BM3DFinal       = core.bm3d.Final
    Crop            = core.std.CropRel
    srclow          = src if srclow is None else srclow
    pad             = padding (src, a+4, a+4, a+4, a+4)
    hflt            = NLMeans (pad, d=0, a=a, s=4, h=h, wref=1.0)
    hif             = Crop (hflt, a+4, a+4, a+4, a+4)
    ref             = BM3DBasic (srclow, sigma=sigma, block_size=block_size, block_step=block_step, group_size=group_size, bm_range=bm_range, bm_step=bm_step)
    lowf            = BM3DFinal (srclow, ref, sigma=sigma, block_size=block_size, block_step=block_step, group_size=group_size, bm_range=bm_range, bm_step=bm_step)
    clip            = hipass (lowf, hif, p=lowpass)
    return clip

def nrfinal (spatial, dif, vmulti, peldif=None, pel=4, tr=6, thsad=400, thscd1=10000, thscd2=255):
    core            = vs.get_core ()
    Expr            = core.std.Expr
    MergeDiff       = core.std.MergeDiff
    MSuper          = core.mvsf.Super
    MDegrainN       = mvmulti.DegrainN
    superclip       = MSuper (dif, pelclip=peldif, rfilter=2, pel=pel, **msuper_args)
    blankd          = Expr ([dif], "0.5")
    coarse          = MDegrainN (blankd, superclip, vmulti, tr=tr, thsad=10000, thscd1=thscd1, thscd2=thscd2, **mdegrain_args)
    fine            = MDegrainN (blankd, superclip, vmulti, tr=tr, thsad=thsad, thscd1=thscd1, thscd2=thscd2, **mdegrain_args)
    newdif          = Expr ([dif, coarse, blankd, fine], ["y z - abs x z - abs > a y ?"])
    clip            = MergeDiff (spatial, newdif)
    return clip

### Sharpeners ###
def delicatesharp (src, relaxed=False, coarse=False, lowpass=16):
    core            = vs.get_core ()
    Median          = core.flt.Median
    Repair          = core.rgsf.Repair
    RemoveGrain     = core.rgsf.RemoveGrain
    Expr            = core.std.Expr
    MakeDiff        = core.std.MakeDiff
    MergeDiff       = core.std.MergeDiff
    def inline (src):
        rg11        = RemoveGrain (src, 11)
        rg11D       = MakeDiff (src, rg11)
        rg11DR      = RemoveGrain (rg11D, 11)
        rg11DD      = Expr ([rg11D, rg11DR], ["x y - x 0.5 - * 0 < 0.5 x y - abs x 0.5 - abs < x y - 0.5 + x ? ?"])
        clip        = MakeDiff (src, rg11DD)
        return clip
    blur            = Repair (Median (src), src, mode=16) if relaxed else Median (src)
    dif             = Expr ([src, blur], ["y x - 3 * 0.5 +"])
    dif             = inline (dif)
    difrg           = RemoveGrain (dif, 11)
    fdif            = Expr ([dif, difrg], ["y 0.5 - abs x 0.5 - abs > y 0.5 ?"])
    merge           = MergeDiff (src, fdif)
    clip            = merge if coarse or relaxed else hipass (src, merge, p=lowpass)
    return clip

def regularsharp (src, median=False, lowpass=16, a=32, h=12.8, thr=0.00390625, elast=None):
    core            = vs.get_core ()
    Median          = core.flt.Median
    resample        = core.fmtc.resample
    MakeDiff        = core.std.MakeDiff
    MergeDiff       = core.std.MergeDiff
    width           = src.width
    height          = src.height
    def inline (clp, n):
        blur        = resample (clp, width*8, height*8, kernel="cubic", a1=1, a2=0, **fmtc_args)
        sharp       = resample (clp, width*8, height*8, kernel="cubic", a1=-1, a2=0, **fmtc_args)
        dif         = MakeDiff (blur, sharp)
        dif         = resample (dif, width*2, height*2, kernel="cubic", a1=-1, a2=0, **fmtc_args)
        clp         = MergeDiff (clp, dif)
        n           = n-1
        return clp if n == 0 else inline (clp, n)
    u2x             = 0 if median else genpelclip (src, pel=2)
    ublur           = 0 if median else inline (u2x, 4)
    cdif            = 0 if median else resample (MakeDiff (u2x, ublur), width, height, -0.5, -0.5, kernel="cubic", a1=-4, a2=-2, **fmtc_args)
    csharp          = 0 if median else halonr (hipass (src, MergeDiff (src, cdif), lowpass), a, h, thr, elast, lowpass//2)
    msharp          = MergeDiff (src, MakeDiff (src, Median (src))) if median else 0
    clip            = msharp if median else csharp
    return clip

def deconvolution (src, loop=2, lowpass=8, a=32, h=12.8, thr=0.00390625, elast=None):
    core            = vs.get_core ()
    Deconv          = core.vcfreq.Sharp
    dcv             = Deconv (src, **deconv_args)
    clip            = hipass (src, dcv, p=lowpass)
    loop            = loop - 1
    return halonr (clip, a, h, thr, elast, lowpass) if loop == 0 else deconvolution (clip, loop, lowpass, a, h, thr, elast)

def sharpfinal (soft, dif, limit, vmulti, peldif=None, pellimit=None, pel=4, tr=6, thsad=400, thscd1=10000, thscd2=255, str=1.00, lowpass=8, a=32, h=6.4, thr=0.00390625, elast=None):
    core            = vs.get_core ()
    NLMeans         = core.knlm.KNLMeansCL
    Repair          = core.rgsf.Repair
    Maximum         = core.flt.Maximum
    Minimum         = core.flt.Minimum
    Expr            = core.std.Expr
    MakeDiff        = core.std.MakeDiff
    MergeDiff       = core.std.MergeDiff
    Crop            = core.std.CropRel
    MSuper          = core.mvsf.Super
    MDegrainN       = mvmulti.DegrainN
    MCompensate     = mvmulti.Compensate
    hintra          = pow (1.2153830693575750222940769140908, str)
    expression      = ["{x} {y} {x} - abs 4 / 1 4 / pow 4 * 1.6 * {str} * {y} {x} - {y} {x} - abs 1.001 + / * + 256 /".format (str=str, x="x 256 *", y="y 256 *")]
    blankd          = Expr ([dif], ["0.5"])
    superdif        = MSuper (dif, pelclip=peldif, rfilter=2, pel=pel, **msuper_args)
    supercmp        = MSuper (limit, pelclip=pellimit, rfilter=2, pel=pel, **msuper_args)
    coarse          = MDegrainN (blankd, superdif, vmulti, tr=tr, thsad=10000, thscd1=thscd1, thscd2=thscd2, **mdegrain_args)
    fine            = MDegrainN (blankd, superdif, vmulti, tr=tr, thsad=thsad, thscd1=thscd1, thscd2=thscd2, **mdegrain_args)
    newdif          = Expr ([dif, coarse, blankd, fine], ["y z - abs x z - abs > a y ?"])
    newdif          = NLMeans (padding (newdif, a+4, a+4, a+4, a+4), d=0, a=a, s=4, h=hintra, wref=1.0, rclip=padding (MergeDiff (soft, newdif), a+4, a+4, a+4, a+4))
    averaged        = MergeDiff (soft, Crop (newdif, a+4, a+4, a+4, a+4))
    comp            = MCompensate (soft, supercmp, vmulti, tr=tr, thsad=thsad, thscd1=thscd1, thscd2=thscd2)
    bright          = Expr ([maxmulti (comp, tr=tr), Maximum (limit)], "x y min")
    dark            = Expr ([minmulti (comp, tr=tr), Minimum (limit)], "x y max")
    clamped         = clamp (averaged, bright, dark, overshoot=0.0, undershoot=0.0)
    amplified       = Expr ([soft, clamped], expression)
    clip            = hipass (soft, Repair (amplified, halonr (amplified, a, h, thr, elast, lowpass), 16), lowpass)
    return clip

### Motion Estimation ###
def getvectors (src, pelclip=None, tr=6, pel=4, thsad=400):
    core            = vs.get_core ()
    MSuper          = core.mvsf.Super
    MAnalyze        = mvmulti.Analyze
    MRecalculate    = mvmulti.Recalculate
    supersoft       = MSuper (src, pelclip=pelclip, rfilter=4, pel=pel, **msuper_args)
    supersharp      = MSuper (src, pelclip=pelclip, rfilter=2, pel=pel, **msuper_args)
    vmulti          = MAnalyze (supersoft, overlap=16, blksize=32, divide=2, tr=tr, dct=7, **manalyze_args)
    vmulti          = MRecalculate (supersoft, vmulti, overlap=8, blksize=16, divide=2, tr=tr, thsad=thsad/2, dct=8, **mrecalculate_args)
    vmulti          = MRecalculate (supersharp, vmulti, overlap=4, blksize=8, divide=2, tr=tr, thsad=thsad/2, dct=6, **mrecalculate_args)
    vmulti          = MRecalculate (supersharp, vmulti, overlap=2, blksize=4, divide=0, tr=tr, thsad=thsad/2, dct=5, **mrecalculate_args)
    return vmulti

def genpelclip (src, src2=None, pel=4):
    core            = vs.get_core ()
    NNEDI           = core.nnedi3.nnedi3
    MakeDiff        = core.std.MakeDiff
    Transpose       = core.std.Transpose
    u2x             = Transpose (NNEDI (Transpose (NNEDI (src, **nnedi_args)), **nnedi_args))
    u2x2            = 0 if src2 is None else Transpose (NNEDI (Transpose (NNEDI (src2, **nnedi_args)), **nnedi_args))
    u4x             = Transpose (NNEDI (Transpose (NNEDI (u2x, **nnedi_args)), **nnedi_args))
    u4x2            = 0 if src2 is None else Transpose (NNEDI (Transpose (NNEDI (u2x2, **nnedi_args)), **nnedi_args))
    dif             = 0 if src2 is None else (MakeDiff (u2x, u2x2) if pel == 2 else MakeDiff (u4x, u4x2))
    clip            = (u2x if pel == 2 else u4x) if src2 is None else dif
    return clip
    
### Resizing ###
def resizenr (src, w=None, h=None, sx=0, sy=0, sw=0, sh=0, kernel="spline", taps=4, a1=None, a2=None, a3=None, center=True):
    core            = vs.get_core ()
    resample        = core.fmtc.resample
    Repair          = core.rgsf.Repair
    BlankClip       = core.std.BlankClip
    MaskedMerge     = core.std.MaskedMerge
    w               = src.width if w is None else w
    h               = src.height if h is None else h
    sr_h            = w / src.width
    sr_v            = h / src.height
    sr_up           = max (sr_h, sr_v)
    sr_dw           = 1.0 / min (sr_h, sr_v)
    sr              = max (sr_up, sr_dw)
    nrb             = (sr > 2.5)
    nrf             = (sr < 2.5 + 1.0)
    nrr             = min (sr - 2.5, 1.0) if nrb else 1.0
    nrv             = (1.0 - nrr) if nrb else 0.0
    nrm             = BlankClip (clip=src, width=w, height=h, color=nrv) if nrb and nrf else 0
    main            = resample (src, w=w, h=h, sx=sx, sy=sy, sw=sw, sh=sh, kernel=kernel, taps=taps, a1=a1, a2=a2, a3=a3, center=center, **fmtc_args)
    nrng            = resample (src, w=w, h=h, sx=sx, sy=sy, sw=sw, sh=sh, kernel="gauss", a1=100, center=center, **fmtc_args) if nrf else main
    clip            = Repair (main, nrng, 1) if nrf else main
    clip            = MaskedMerge (main, clip, nrm) if nrf and nrb else clip
    return clip

### Chroma Reconstruction ###
def subtofull (src, cplace="mpeg2"):
    core            = vs.get_core ()
    resample        = core.fmtc.resample
    Expr            = core.std.Expr
    ShufflePlanes   = core.std.ShufflePlanes
    factor_w        = src.format.subsampling_w
    factor_h        = src.format.subsampling_h
    def inline_pel_422 (src):
        NNEDI       = core.nnedi3.nnedi3
        Transpose   = core.std.Transpose
        clip        = Transpose (NNEDI (Transpose (src), **nnedi_args))
        return clip
    srcy            = ShufflePlanes (src, 0, vs.GRAY)
    srcu            = Expr (ShufflePlanes (src, 1, vs.GRAY), "x 0.5 +")
    srcv            = Expr (ShufflePlanes (src, 2, vs.GRAY), "x 0.5 +")
    if factor_w == 1 and factor_h == 1:
       unew         = genpelclip (srcu, pel=2)
       vnew         = genpelclip (srcv, pel=2)
       if cplace == "mpeg2":
          ushift    = resample (unew, sx=0, sy=-0.5, kernel="spline", taps=6, **fmtc_args)
          vshift    = resample (vnew, sx=0, sy=-0.5, kernel="spline", taps=6, **fmtc_args)
       elif cplace == "mpeg1":
          ushift    = resample (unew, sx=-0.5, sy=-0.5, kernel="spline", taps=6, **fmtc_args)
          vshift    = resample (vnew, sx=-0.5, sy=-0.5, kernel="spline", taps=6, **fmtc_args)
       else:
          raise ValueError ("wtf?")
    if factor_w == 1 and factor_h == 0:
       unew         = inline_pel_422 (srcu)
       vnew         = inline_pel_422 (srcv)
       ushift       = resample (unew, sx=1, sy=0, kernel="spline", taps=6, **fmtc_args)
       vshift       = resample (vnew, sx=1, sy=0, kernel="spline", taps=6, **fmtc_args)
    ufinal          = Expr (ushift, "x 0.5 -")
    vfinal          = Expr (vshift, "x 0.5 -")
    clip            = ShufflePlanes ([srcy, ufinal, vfinal], [0, 0, 0], vs.YUV)
    return clip

def fulltonative (src, a=32, h=6.4, lowpass=6, mode=0):
    core            = vs.get_core ()
    NLMeans         = core.knlm.KNLMeansCL
    resample        = core.fmtc.resample
    Expr            = core.std.Expr
    MakeDiff        = core.std.MakeDiff
    MergeDiff       = core.std.MergeDiff
    ShufflePlanes   = core.std.ShufflePlanes
    Crop            = core.std.CropRel
    def gauss_h (src, p=30):
        upsmp       = resample (src, src.width * 2, src.height, kernel="gauss", a1=100, **fmtc_args)
        clip        = resample (upsmp, src.width, src.height, kernel="gauss", a1=p, **fmtc_args)
        return clip
    pad             = padding (src, a, a, a, a)
    srcy            = ShufflePlanes (pad, 0, vs.GRAY)
    srcu            = ShufflePlanes (pad, 1, vs.GRAY)
    srcv            = ShufflePlanes (pad, 2, vs.GRAY)
    u4x             = genpelclip (pad, pel=4)
    luma            = ShufflePlanes (u4x, 0, vs.GRAY)
    u               = Expr (ShufflePlanes (u4x, 1, vs.GRAY), "x 0.5 +")
    v               = Expr (ShufflePlanes (u4x, 2, vs.GRAY), "x 0.5 +")
    unew            = resizenr (NLMeans (u, d=0, a=a, s=0, h=h, wref=1.0, rclip=luma), pad.width, pad.height, sx=-1.5, sy=-1.5, kernel="spline", taps=6)
    vnew            = resizenr (NLMeans (v, d=0, a=a, s=0, h=h, wref=1.0, rclip=luma), pad.width, pad.height, sx=-1.5, sy=-1.5, kernel="spline", taps=6)
    uhi             = MakeDiff (unew, gauss_h (unew, p=lowpass)) if mode else MakeDiff (unew, gauss (unew, p=lowpass))
    vhi             = MakeDiff (vnew, gauss_h (vnew, p=lowpass)) if mode else MakeDiff (vnew, gauss (vnew, p=lowpass))
    ufinal          = Expr (MergeDiff (Expr (srcu, "x 0.5 +"), uhi), "x 0.5 -")
    vfinal          = Expr (MergeDiff (Expr (srcv, "x 0.5 +"), vhi), "x 0.5 -")
    merge           = ShufflePlanes ([srcy, ufinal, vfinal], [0, 0, 0], vs.YUV)
    clip            = Crop (merge, a, a, a, a)
    return clip

### Sigmoidal ###
def build_sigmoid_expr (string, inv=False, thr=0.5, cont=6.5):
    x1m0            = "1 {thr} 1 - {cont} * exp 1 + / 1 {cont} {thr} * exp 1 + / -".format (thr=thr, cont=cont)
    x0              = "1 {cont} {thr} * exp 1 + /".format (thr=thr, cont=cont)
    if inv:
       expr         = "{thr} 1 " + string + " {x1m0} * {x0} + 0.000001 max / 1 - 0.000001 max log {cont} / -".format (x1m0=x1m0, x0=x0, thr=thr, cont=cont)
    else:
       expr         = "1 1 {cont} {thr} " + string + " - * exp + / {x0} - {x1m0} /".format (x1m0=x1m0, x0=x0, thr=thr, cont=cont)
    return expr.format (thr=thr, cont=cont)

def sigmoid_direct (src, thr=0.5, cont=6.5):
    core            = vs.get_core ()
    Expr            = core.std.Expr
    clip            = Expr (src, build_sigmoid_expr ("x", False, thr, cont))
    return clip

def sigmoid_inverse (src, thr=0.5, cont=6.5):
    core            = vs.get_core ()
    Expr            = core.std.Expr
    clip            = Expr (src, build_sigmoid_expr ("x", True, thr, cont))
    return clip
