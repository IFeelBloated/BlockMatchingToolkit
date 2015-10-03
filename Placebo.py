import vapoursynth as vs
import math
import mvmulti

### Shared ###
def gauss (src, p=30):
    core   = vs.get_core ()
    upsmp  = core.fmtc.resample (src, src.width * 2, src.height * 2, kernel="gauss", a1=100, fulls=True, fulld=True)
    clip   = core.fmtc.resample (upsmp, src.width, src.height, kernel="gauss", a1=p, fulls=True, fulld=True)
    return clip

def padding (src, left=0, right=0, top=0, bottom=0):
    core   = vs.get_core ()
    w      = src.width
    h      = src.height
    clip   = core.fmtc.resample (src, w+left+right, h+top+bottom, -left, -top, w+left+right, h+top+bottom, kernel="point", fulls=True, fulld=True)
    return clip

def xymax (src1, src2):
    core   = vs.get_core ()
    clip   = core.std.Expr ([src1, src2], ["x y > x y ?"])
    return clip

def xymin (src1, src2):
    core   = vs.get_core ()
    clip   = core.std.Expr ([src1, src2], ["x y > y x ?"])
    return clip

def max_dif (src1, src2, ref):
    core = vs.get_core ()
    clip = core.std.Expr ([src1, src2, ref], ["x z - abs y z - abs > x y ?"])
    return clip

def min_dif (src1, src2, ref):
    core = vs.get_core ()
    clip = core.std.Expr ([src1, src2, ref], ["x z - abs y z - abs > y x ?"])
    return clip

def thr_merge (flt, src, ref=None, thr=0.0009765625, elast=None):
    core  = vs.get_core ()
    ref   = src if ref is None else ref
    elast = thr / 2 if elast is None else elast
    BExp  = ["x {thr} {elast} + z - 2 {elast} * / * y {elast} z + {thr} - 2 {elast} * / * +".format (thr=thr, elast=elast)]
    BDif  = core.std.Expr (src, "0.0")
    PDif  = core.std.Expr ([flt, src], "x y - 0.0 max")
    PRef  = core.std.Expr ([flt, ref], "x y - 0.0 max")
    PBLD  = core.std.Expr ([PDif, BDif, PRef], BExp)
    NDif  = core.std.Expr ([flt, src], "y x - 0.0 max")
    NRef  = core.std.Expr ([flt, ref], "y x - 0.0 max")
    NBLD  = core.std.Expr ([NDif, BDif, NRef], BExp)
    BLDD  = core.std.Expr ([PBLD, NBLD], ["x y - 0.5 +"])
    BLD   = core.std.MergeDiff (src, BLDD)
    UDN   = core.std.Expr ([flt, ref, BLD], ["x y - abs {thr} {elast} - > z x ?".format (thr=thr, elast=elast)])
    clip  = core.std.Expr ([flt, ref, UDN, src], ["x y - abs {thr} {elast} + < z a ?".format (thr=thr, elast=elast)])
    return clip

### Denoise ###
def halonr (src, a=32, h=6.4, thr=0.00390625, elast=None):
    core    = vs.get_core ()
    elast   = thr / 6 if elast is None else elast
    pad     = padding (src, a, a, a, a)
    Clean   = core.knlm.KNLMeansCL (pad, d=0, a=a, s=0, h=h)
    EMask   = core.tcanny.TCanny (Clean, sigma=1.5, mode=1, op=0)
    EMask   = core.std.Expr (EMask, "x 0.24 - 3.2 * 0.0 max 1.0 min")
    EMask   = core.flt.Maximum (EMask)
    EMask   = core.flt.Inflate (EMask)
    MRG     = core.std.MaskedMerge (pad, Clean, EMask)
    Final   = thr_merge (pad, MRG, thr=thr, elast=elast)
    clip    = core.std.CropRel (Final, a, a, a, a)
    return clip

def ringnr (src, a=32, h=6.4, divide=4, thr=0.03125, elast=None, lowpass=8):
    core    = vs.get_core ()
    elast   = thr / 8 if elast is None else elast
    hfine   = h / divide
    hcoarse = h / (divide - 1)
    pad     = padding (src, a+1, a+1, a+1, a+1)
    def inline (src, h, n):
        flt = core.knlm.KNLMeansCL (src, d=0, a=a, s=1, h=h)
        n   = n - 1
        return flt if n == 0 else inline (flt, h, n)
    coarse  = core.knlm.KNLMeansCL (pad, d=0, a=a, s=1, h=h)
    fine    = inline (pad, hfine, divide)
    fine    = min_dif (fine, coarse, pad)
    coarse  = core.rgsf.Repair (fine, coarse, mode=1)
    ref     = thr_merge (fine, coarse, thr=thr, elast=elast)
    clean   = core.knlm.KNLMeansCL (pad, d=0, a=a, s=1, h=hcoarse, rclip=ref)
    dif     = core.std.MakeDiff (pad, clean)
    dif     = core.knlm.KNLMeansCL (dif, d=0, a=a, s=1, h=hcoarse, rclip=clean)
    fine    = core.std.MergeDiff (clean, dif)
    coarse  = core.rgsf.Repair (fine, clean, mode=9)
    hif     = thr_merge (fine, coarse, thr=thr, elast=elast)
    mrg     = core.std.MergeDiff (gauss (pad, p=lowpass), core.std.MakeDiff (hif, gauss (hif, p=lowpass)))
    clip    = core.std.CropRel (mrg, a+1, a+1, a+1, a+1)
    return clip

def spatialnr (src, a=32, h=1.2, sbsize=33, sstring="0.0:16.0 0.48:8.0 0.64:0.5 1.0:0.0", lowpass=12):
    core    = vs.get_core ()
    hfloor  = math.floor (h)
    hcoarse = (hfloor * 2.1445390776269709272958045814962) + (h - hfloor)
    pad     = padding (src, a+4, a+4, a+4, a+4)
    NLM     = core.knlm.KNLMeansCL (pad, d=0, a=a, s=4, h=h)
    DFT     = core.dfttest.DFTTest (pad, sstring=sstring, sbsize=sbsize, sosize=0, smode=0, tosize=0, tbsize=1, tmode=0)
    MRG     = core.std.MakeDiff (NLM, gauss (NLM, p=lowpass)).std.MergeDiff (gauss (DFT, p=lowpass))
    Dif     = core.std.MakeDiff (pad, MRG)
    Dif     = core.knlm.KNLMeansCL (Dif, d=0, a=a, s=4, h=hcoarse, rclip=MRG)
    Final   = core.std.MergeDiff (MRG, Dif)
    clip    = core.std.CropRel (Final, a+4, a+4, a+4, a+4)
    return clip

def nrfinal (spatial, dif, peldif, vmulti, pel=4, tr=6, thsad=10000, thscd1=10000, thscd2=255, repmode=13):
    core      = vs.get_core ()
    superclip = core.mvsf.Super (dif, pel=pel, chroma=False, hpad=32, vpad=32, pelclip=peldif, sharp=2, rfilter=2, levels=0)
    blankd    = core.std.Expr ([dif], "0.5")
    mc        = mvmulti.DegrainN (blankd, superclip, vmulti, tr=tr, thsad=thsad, thscd1=thscd1, thscd2=thscd2, plane=0, limit=1.0)
    NR        = core.std.MergeDiff (spatial, mc)
    repclip   = core.std.MergeDiff (spatial, dif)
    clip      = core.rgsf.Repair (NR, repclip, repmode)
    return clip

### ME & MC ###
def genpelclip (src, pel=4):
    core   = vs.get_core ()
    w      = src.width
    h      = src.height
    if pel == 2:
       clip = ediresample (src, w*2, h*2, sx=0.25, sy=0.25, noring=True, kernel_u="spline", kernel_d="spline", taps=6, fulls=True, fulld=True)
    elif pel == 4:
       clip = ediresample (src, w*4, h*4, sx=0.375, sy=0.375, noring=True, kernel_u="spline", kernel_d="spline", taps=6, fulls=True, fulld=True)
    else:
       clip = 0
    return clip

def getvectors (src, pelclip, tr=6, pel=4, dct=5, thsad=400):
    core       = vs.get_core ()
    supersoft  = core.mvsf.Super (src, pel=pel, chroma=False, hpad=32, vpad=32, pelclip=pelclip, sharp=2, rfilter=4, levels=0)
    supersharp = core.mvsf.Super (src, pel=pel, chroma=False, hpad=32, vpad=32, pelclip=pelclip, sharp=2, rfilter=2, levels=0)
    vmulti     = mvmulti.Analyze (supersoft, overlap=16, blksize=32, search=3, chroma=False, truemotion=True, tr=tr, trymany=True, searchparam=16, pelsearch=16, dct=dct, levels=0, divide=2, badrange=-24)
    vmulti     = mvmulti.Recalculate (supersoft, vmulti, tr=tr, overlap=8, blksize=16, thsad=thsad/2, chroma=False, truemotion=True, search=3, searchparam=16, dct=dct, smooth=1, divide=2)
    vmulti     = mvmulti.Recalculate (supersharp, vmulti, tr=tr, overlap=4, blksize=8, thsad=thsad/2, chroma=False, truemotion=True, search=3, searchparam=16, dct=dct, smooth=1, divide=2)
    vmulti     = mvmulti.Recalculate (supersharp, vmulti, tr=tr, overlap=2, blksize=4, thsad=thsad/2, chroma=False, truemotion=True, search=3, searchparam=16, dct=dct, smooth=1, divide=0)
    return vmulti

### Resize ###
def resizenr (src, w=None, h=None, sx=0, sy=0, sw=0, sh=0, kernel="spline36", kernelh=None, kernelv=None, fh=1, fv=1, taps=4, a1=None, a2=None, a3=None, kovrspl=1, cnorm=True, center=True, fulls=None, fulld=None, cplace="mpeg2", invks=False, invkstaps=4, noring=True):
    core    = vs.get_core ()
    scss    = GetCSS (src)
    Gray    = scss == "GRAY"
    w       = src.width if w is None else w
    h       = src.height if h is None else h
    kernelh = kernel if kernelh is None else kernelh
    kernelv = kernel if kernelv is None else kernelv
    sr_h    = float (w / src.width)
    sr_v    = float (h / src.height)
    sr_up   = max (sr_h, sr_v)
    sr_dw   = 1.0 / min (sr_h, sr_v)
    sr      = max (sr_up, sr_dw)
    thr   = 2.5
    nrb   = (sr > thr)
    nrf   = (sr < thr + 1.0 and noring)
    nrr   = min (sr - thr, 1.0) if nrb else 1.0
    if Gray:
       nrv = 1.0 - nrr if nrb else 0
    else:
       nrv = [1.0 - nrr, 1.0 - nrr, 1.0 - nrr] if nrb else [0, 0, 0]
    nrm   = core.std.BlankClip (clip=src, width=w, height=h, color=nrv) if nrb and nrf else 0
    main  = core.fmtc.resample (src, w=w, h=h, sx=sx, sy=sy, sw=sw, sh=sh, kernel=kernel, kernelh=kernelh, kernelv=kernelv, fh=fh, fv=fv, taps=taps, a1=a1, a2=a2, a3=a3, kovrspl=kovrspl, cnorm=cnorm, center=center, fulls=fulls, fulld=fulld, cplace=cplace, invks=invks, invkstaps=invkstaps)
    nrng  = core.fmtc.resample (src, w=w, h=h, sx=sx, sy=sy, sw=sw, sh=sh, kernel="gauss", a1=100, center=center, fulls=fulls, fulld=fulld, cplace=cplace) if nrf else main
    clip  = core.rgsf.Repair (main, nrng, 1) if nrf else main
    clip  = core.std.MaskedMerge (main, clip, nrm) if nrf and nrb else clip
    return clip

def GetCSS (src):
    core   = vs.get_core ()
    if src.format.id == vs.GRAYS:
       css = "GRAY"
    else:
       Y   = core.std.ShufflePlanes (src, planes=0, colorfamily=vs.GRAY)
       C   = core.std.ShufflePlanes (src, planes=1, colorfamily=vs.GRAY)
       YW  = Y.width
       YH  = Y.height
       CW  = C.width
       CH  = C.height
       if CH == int (0.5 * YH):
          css = "420"
       elif CW == int (0.5 * YW):
          css = "422"
       else:
          css = "444"
    return css

def ediresample (src, w=None, h=None, sx=0, sy=0, sw=None, sh=None, kernel_u="spline64", kernel_d="bicubic", taps=4, a1=-1, a2=0, a3=None, css=None, fulls=False, fulld=None, cplace="mpeg2", nsize=0, nns=4, qual=2, etype=0, pscrn=1, noring=False, ratiothr=1.0):
    core   = vs.get_core ()
    scss   = GetCSS (src)
    css    = scss if css is None else css
    fulld  = fulls if fulld is None else fulld
    Gray   = scss == "GRAY" or css == "GRAY"
    scss   = "444" if Gray else scss
    ocss   = "444" if Gray else css
    InterK = resizenr if noring else core.fmtc.resample
    Y      = 3
    U      = 1 if Gray else 3
    V      = 1 if Gray else 3
    Ut     = U == 3
    Vt     = V == 3
    cw     = src.width
    ch     = src.height
    cwc    = cw if scss == "444" else cw // 2
    chc    = ch // 2 if scss == "420" else ch
    ow     = cw if w is None else w
    oh     = ch if h is None else h
    owc    = ow if ocss == "444" else ow // 2
    ohc    = oh // 2 if ocss == "420" else oh
    sw     = cw if sw is None else sw
    sh     = ch if sh is None else sh
    prel   = int (sx / 2) * 2
    pret   = int (sy / 2) * 2
    prer   = int ((-cw + sx + sw if sw > 0 else sw) / 2) * 2
    preb   = int ((-ch + sy + sh if sh > 0 else sh) / 2) * 2
    prew   = cw - prel + prer
    preh   = ch - pret + preb
    if scss == "444":
       cwmod2 = int (cw / 2) * 2 == cw
       pwmod2 = int (prew / 2) * 2 == prew
       wpre   = prew < cw
       prel   = prel if wpre else 0
       prer   = (prer if pwmod2 else prer + 1) if wpre else (0 if cwmod2 else 1)
       prew   = cw - prel + prer
       wpre   = prew < cw or cwmod2 == False
    else:
       cwmod4 = int (cw / 4) * 4 == cw
       pwmod4 = int (prew / 4) * 4 == prew
       wpre   = prew < cw
       prel   = prel if wpre else 0
       prer   = (prer if pwmod4 else prer + 2) if wpre else (0 if cwmod4 else 2)
       prew   = cw - prel + prer
       wpre   = prew < cw or cwmod4 == False
    if scss == "420":
       chmod4 = int (ch / 4) * 4 == ch
       phmod4 = int (preh / 4) * 4 == preh
       hpre   = preh < ch
       pret   = pret if hpre else 0
       preb   = (preb if phmod4 else preb + 2) if hpre else (0 if chmod4 else 2)
       preh   = ch - pret + preb
       hpre   = preh < ch or chmod4 == False
    else:
       chmod2 = int (ch / 2) * 2 == ch
       phmod2 = int (preh / 2) * 2 == preh
       hpre   = preh < ch
       pret   = pret if hpre else 0
       preb   = (preb if phmod2 else preb + 1) if hpre else (0 if chmod2 else 1)
       preh   = ch - pret + preb
       hpre   = preh < ch or chmod2 == False
    sw     = cw - sx + sw if sw <= 0 else sw
    sh     = ch - sy + sh if sh <= 0 else sh
    sx     = sx - prel if wpre else sx
    sy     = sy - pret if hpre else sy
    sxc    = sx if scss == "444" else float (sx / 2)
    syc    = float (sy / 2) if scss == "420" else sy
    swc    = sw if scss == "444" else float (sw / 2)
    shc    = float (sh / 2) if scss == "420" else sh
    yhratio = float (ow / sw)
    yvratio = float (oh / sh)
    chratio = float (owc / swc)
    cvratio = float (ohc / shc)
    enable   = yhratio != 1 or yvratio != 1 or chratio != 1 or cvratio != 1 or sw != int (sw) or sh != int (sh) or swc != int (swc) or shc != int (shc) or sx != int (sx) or sy != int (sy) or sxc != int (sxc) or syc != int (syc)
    yhct     = math.ceil ((math.log (yhratio / ratiothr)) / math.log (2)) if yhratio > ratiothr else 0
    yhrf     = int (math.pow (2, yhct))
    yrhratio = yhratio / yhrf
    yvct     = math.ceil ((math.log (yvratio / ratiothr)) / math.log (2)) if yvratio > ratiothr else 0
    yvrf     = int (math.pow (2, yvct))
    yrvratio = yvratio / yvrf
    chct     = math.ceil ((math.log (chratio / ratiothr)) / math.log (2)) if chratio > ratiothr else 0
    chrf     = int (math.pow (2, chct))
    crhratio = chratio / chrf
    cvct     = math.ceil ((math.log (cvratio / ratiothr)) / math.log (2)) if cvratio > ratiothr else 0
    cvrf     = int (math.pow (2, cvct))
    crvratio = cvratio / cvrf
    noediy = yhct <= 0 and yvct <= 0
    noedic = chct <= 0 and cvct <= 0
    noedi  = noediy or noedic
    Yedit  = noediy == False
    Uedit  = Ut and noedic == False
    Vedit  = Vt and noedic == False
    edit   = Yedit or Uedit or Vedit
    mixed  = False if edit == False or enable == False else True
    yhchift = 0.5 if yhrf >= 2 else 0
    yvchift = 0.5 if yvrf >= 2 else 0
    yhfix   = -yhchift
    yvfix   = -yvchift
    chshift = ((0.5 if chrf >= 2 else 0) if scss == "444" else ((0.5 if chrf >= 2 else 0) if cplace == "mpeg1" else (0.5 - float (chrf / 4) if chrf >= 2 else -0.25))) if ocss == "444" else (((0.5 if chrf >= 2 else 0) if cplace == "mpeg1" else (0.75 if chrf >= 2 else 0.25)) if scss == "444" else ((0.5 if chrf >= 2 else 0) if cplace == "mpeg1" else (0.75 - float (chrf / 4) if chrf >= 2 else 0)))
    cvshift = 0.5 if cvrf >= 2 else 0
    chfix   = -chshift
    cvfix   = -cvshift
    cphfixe = 0 if ocss == "444" else (0 if cplace == "mpeg1" else 0.25 - float (0.25 / crhratio))
    cphfix  = (0 if scss == "444" else (0 if cplace == "mpeg1" else 0.25)) if ocss=="444" else ((0 if cplace == "mpeg1" else -0.5) if scss == "444" else (0 if cplace == "mpeg1" else 0.25 - float (0.25 / chratio)))
    Luma    = core.std.ShufflePlanes (src, planes=0, colorfamily=vs.GRAY) if GetCSS (src) != "GRAY" else src
    input   = src if Gray == False else Luma
    input   = core.fmtc.resample (input, (prew if wpre else cw), (preh if hpre else ch), (prel if wpre else 0), (pret if hpre else 0), (prew if wpre else cw), (preh if hpre else ch), kernel="point", fulls=fulls, fulld=fulls) if wpre or hpre else input
    if enable == False and edit == False:
       0
    elif yhct == chct and yvct == cvct and scss == "420":
       edgeedi  = EDInter (input, yvct, yhct, 1, 1, [0, 1, 2], nsize, nns, qual, etype, pscrn)
       edgeediY = core.std.ShufflePlanes (edgeedi, planes=0, colorfamily=vs.GRAY)
       edgeediU = core.std.Expr (core.std.ShufflePlanes (edgeedi, planes=1, colorfamily=vs.GRAY), "x 0.5 +")
       edgeediV = core.std.Expr (core.std.ShufflePlanes (edgeedi, planes=2, colorfamily=vs.GRAY), "x 0.5 +")
    else:
       edgeediY = core.std.ShufflePlanes (input, planes=0, colorfamily=vs.GRAY) if GetCSS (src) != "GRAY" else input
       edgeediY = edgeediY if Yedit == False else EDInter (edgeediY, yvct, yhct, 1, 1, 0, nsize, nns, qual, etype, pscrn)
       edgeediU = core.std.Expr (core.std.ShufflePlanes (input, planes=1, colorfamily=vs.GRAY), "x 0.5 +") if Gray == False else 0
       edgeediU = 0 if Uedit == False else EDInter (edgeediU, cvct, chct, 1, 1, 0, nsize, nns, qual, etype, pscrn)
       edgeediV = core.std.Expr (core.std.ShufflePlanes (input, planes=2, colorfamily=vs.GRAY), "x 0.5 +") if Gray == False else 0
       edgeediV = 0 if Vedit == False else EDInter (edgeediV, cvct, chct, 1, 1, 0, nsize, nns, qual, etype, pscrn)
    yrh = yrhratio > ratiothr
    yrv = yrvratio > ratiothr
    crh = crhratio > ratiothr
    crv = crvratio > ratiothr
    if enable == False and edit == False:
       0
    else:
       edgeY   = edgeediY if Yedit == False else core.fmtc.resample (edgeediY, ow, oh, (sx * yhrf + yhfix), (sy * yvrf + yvfix), (sw * yhrf), (sh * yvrf), kernelh=(kernel_u if yrh else kernel_d), kernelv=(kernel_u if yrv else kernel_d), taps=taps, a1=a1, a2=a2, a3=a3, fulls=fulls, fulld=fulls)
       edgeU   = 0 if Uedit == False else core.std.Expr (core.fmtc.resample (edgeediU, owc, ohc, (sxc * chrf + chfix + cphfixe), (syc * cvrf + cvfix), (swc * chrf), (shc * cvrf), kernelh=(kernel_u if crh else kernel_d), kernelv=(kernel_u if crv else kernel_d), taps=taps, a1=a1, a2=a2, a3=a3, fulls=fulls, fulld=fulls), "x 0.5 -")	
       edgeV   = 0 if Vedit == False else core.std.Expr (core.fmtc.resample (edgeediV, owc, ohc, (sxc * chrf + chfix + cphfixe), (syc * cvrf + cvfix), (swc * chrf), (shc * cvrf), kernelh=(kernel_u if crh else kernel_d), kernelv=(kernel_u if crv else kernel_d), taps=taps, a1=a1, a2=a2, a3=a3, fulls=fulls, fulld=fulls), "x 0.5 -")
       edge    = edgeY if Gray or Uedit == False or Vedit == False else core.std.ShufflePlanes ([edgeY, edgeU, edgeV], planes=[0, 0, 0], colorfamily=vs.YUV)
    yh = yhratio > ratiothr
    yv = yvratio > ratiothr
    ch = chratio > ratiothr
    cv = cvratio > ratiothr
    if enable == False and (mixed == False or (Yedit and Uedit and Vedit)):
       0
    elif yhratio == chratio and yvratio == cvratio and (mixed == False or (Yedit and Uedit and Vedit)):
       flat = InterK (input, ow, oh, sx, sy, sw, sh, cplace=cplace, kernelh=(kernel_u if yh else kernel_d), kernelv=(kernel_u if yv else kernel_d), taps=taps, a1=a1, a2=a2, a3=a3, fulls=fulls, fulld=fulls)
    else:
       flatY = core.std.ShufflePlanes (input, planes=0, colorfamily=vs.GRAY) if GetCSS (src) != "GRAY" else input
       flatY = InterK (flatY, ow, oh, sx, sy, sw, sh, kernelh=(kernel_u if yh else kernel_d), kernelv=(kernel_u if yv else kernel_d), taps=taps, a1=a1, a2=a2, a3=a3, fulls=fulls, fulld=fulls) if mixed or Yedit == False else flatY
       flatU = core.std.Expr (core.std.ShufflePlanes (input, planes=1, colorfamily=vs.GRAY), "x 0.5 +") if (mixed or Uedit == False) and Ut else 0
       flatU = core.std.Expr (InterK (flatU, owc, ohc, (sxc + cphfix), syc, swc, shc, kernelh=(kernel_u if ch else kernel_d), kernelv=(kernel_u if cv else kernel_d), taps=taps, a1=a1, a2=a2, a3=a3, fulls=fulls, fulld=fulls), "x 0.5 -") if (mixed or Uedit == False) and Ut else 0
       flatV = core.std.Expr (core.std.ShufflePlanes (input, planes=2, colorfamily=vs.GRAY), "x 0.5 +") if (mixed or Vedit == False) and Vt else 0
       flatV = core.std.Expr (InterK (flatV, owc, ohc, (sxc + cphfix), syc, swc, shc, kernelh=(kernel_u if ch else kernel_d), kernelv=(kernel_u if cv else kernel_d), taps=taps, a1=a1, a2=a2, a3=a3, fulls=fulls, fulld=fulls), "x 0.5 -") if (mixed or Vedit == False) and Vt else 0	
       flat  = flatY if Gray else core.std.ShufflePlanes ([flatY, flatU, flatV], planes=[0, 0, 0], colorfamily=vs.YUV)
    merge  = edge if edit else flat
    mergeY = core.std.ShufflePlanes (merge, planes=0, colorfamily=vs.GRAY)
    mergeU = core.std.ShufflePlanes (merge, planes=1, colorfamily=vs.GRAY) if Gray == False else 0
    mergeV = core.std.ShufflePlanes (merge, planes=2, colorfamily=vs.GRAY) if Gray == False else 0
    mergeF = core.std.ShufflePlanes ([mergeY, mergeU, mergeV], planes=[0, 0, 0], colorfamily=vs.YUV) if css != "GRAY" else mergeY
    Range  = core.fmtc.bitdepth (mergeF, fulls=fulls, fulld=fulld, dmode=1) if fulls != fulld else mergeF
    return Range

def EDInter (src, vct=1, hct=1, vfield=1, hfield=1, planes=0, nsize=0, nns=4, qual=2, etype=0, pscrn=1, honly=False):
    core = vs.get_core ()
    clip = src
    if hct >= 1:
       clip   = clip if honly else core.std.Transpose (clip)
       clip   = core.nnedi3.nnedi3 (clip, field=hfield, dh=True, planes=planes, nsize=nsize, nns=nns, qual=qual, etype=etype, pscrn=pscrn)
       hct    = hct - 1
       honly  = hct >= 1
       hfield = 0
       clip   = clip if honly else core.std.Transpose (clip)
    else:
       0
    if vct >= 1 and honly == False:
       clip   = core.nnedi3.nnedi3(clip, field=vfield, dh=True, planes=planes, nsize=nsize, nns=nns, qual=qual, etype=etype, pscrn=pscrn)
       vct    = vct - 1
       vfield = 0
    else: 
       0
    clip = clip if vct <= 0 and hct <= 0 else EDInter (clip, vct, hct, vfield, hfield, planes, nsize, nns, qual, etype, pscrn, honly)
    return clip

### Curves ###
def build_sigmoid_expr (string, inv=False, thr=0.5, cont=6.5):
    x1m0 = "1 {thr} 1 - {cont} * exp 1 + / 1 {cont} {thr} * exp 1 + / -".format (thr=thr, cont=cont)
    x0   = "1 {cont} {thr} * exp 1 + /".format (thr=thr, cont=cont)
    if inv:
       expr = "{thr} 1 " + string + " {x1m0} * {x0} + 0.000001 max / 1 - 0.000001 max log {cont} / -".format (x1m0=x1m0, x0=x0, thr=thr, cont=cont)
    else:
       expr = "1 1 {cont} {thr} " + string + " - * exp + / {x0} - {x1m0} /".format (x1m0=x1m0, x0=x0, thr=thr, cont=cont)
    return expr.format (thr=thr, cont=cont)

def sigmoid_direct (src, thr=0.5, cont=6.5):
    core = vs.get_core ()
    expr = build_sigmoid_expr ("x", False, thr, cont)
    clip = core.std.Expr ([src], [expr])
    return clip

def sigmoid_inverse (src, thr=0.5, cont=6.5):
    core = vs.get_core ()
    expr = build_sigmoid_expr ("x", True, thr, cont)
    clip = core.std.Expr ([src], [expr])
    return clip
