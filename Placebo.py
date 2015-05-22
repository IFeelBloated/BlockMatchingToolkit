import vapoursynth as vs
import math

def sbr16 (src):
    core   = vs.get_core ()
    rg11   = core.rgvs.RemoveGrain (src, 11)
    rg11D  = core.std.MakeDiff (src, rg11)
    rg11DR = core.rgvs.RemoveGrain (rg11D, 11)
    rg11DD = core.std.Expr ([rg11D, rg11DR], ["x y - x 32768 - * 0 < 32768 x y - abs x 32768 - abs < x y - 32768 + x ? ?"])
    clip   = core.std.MakeDiff (src, rg11DD)
    return clip

def shrinksharp16 (src, str=1.0):
    core   = vs.get_core ()
    blur   = core.ctmf.CTMF (src, radius=1)
    blur   = core.rgvs.Repair (blur, src, 16)
    dif    = core.std.Expr ([src, blur], ["y x - {str} * 3 * 32768 +".format (str=str)])
    dif    = sbr16 (dif)
    difrg  = core.rgvs.RemoveGrain (dif, 11)
    fdif   = core.std.Expr ([dif, difrg], ["y 32768 - abs x 32768 - abs > y 32768 ?"])
    clip   = core.std.MergeDiff (src, fdif)
    return clip

def sharplimit16 (src, sharp, str=1.0):
    core   = vs.get_core ()
    clip   = core.std.Expr ([src, sharp], ["x 256 / y 256 / x 256 / - abs 4 / log 1 4 / * exp 4 * {str} * y 256 / x 256 / - y 256 / x 256 / - abs 1.001 + / * + 256 *".format (str=str)])
    return clip

def gauss16 (src, p=30):
    core   = vs.get_core ()
    float  = core.fmtc.bitdepth (src, bits=32, flt=True, fulls=True, fulld=True, dmode=1)
    upsmp  = core.fmtc.resample (float, float.width * 2, float.height * 2, kernel="gauss", a1=100, fulls=True, fulld=True)
    gauss  = core.fmtc.resample (upsmp, float.width, float.height, kernel="gauss", a1=p, fulls=True, fulld=True)
    clip   = core.fmtc.bitdepth (gauss, bits=16, flt=False, fulls=True, fulld=True, dmode=1)
    return clip

def xymax16 (src1, src2):
    core   = vs.get_core ()
    clip   = core.std.Expr ([src1, src2], ["x y > x y ?"])
    return clip

def xymin16 (src1, src2):
    core   = vs.get_core ()
    clip   = core.std.Expr ([src1, src2], ["x y > y x ?"])
    return clip

def maxmulti16 (src, start=None, a=2, tr=6):
    core   = vs.get_core ()
    start  = xymax16 (core.std.SelectEvery (src, tr*2+1, 0), core.std.SelectEvery (src, tr*2+1, 1)) if start is None else start
    max    = xymax16 (start, core.std.SelectEvery (src, tr*2+1, a))
    a      = a+1
    clip   = max if a == tr*2+1 else maxmulti16 (src, start=max, tr=tr, a=a)
    return clip

def minmulti16 (src, start=None, a=2, tr=6):
    core   = vs.get_core ()
    start  = xymin16 (core.std.SelectEvery (src, tr*2+1, 0), core.std.SelectEvery (src, tr*2+1, 1)) if start is None else start
    min    = xymin16 (start, core.std.SelectEvery (src, tr*2+1, a))
    a      = a+1
    clip   = min if a == tr*2+1 else minmulti16 (src, start=min, tr=tr, a=a)
    return clip

def spatialnr16 (src, NLM, DFT, lowpass=16, safe=1):
    core   = vs.get_core ()
    DFTL   = gauss16 (DFT, p=lowpass)
    NLML   = gauss16 (NLM, p=lowpass)
    gauss  = gauss16 (src, p=safe)
    NLMH   = core.std.MakeDiff (NLM, NLML)
    HNR    = core.std.MergeDiff (DFTL, NLMH)
    HNRL   = gauss16 (HNR, p=safe)
    HNRH   = core.std.MakeDiff (HNR, HNRL)
    clip   = core.std.MergeDiff (gauss, HNRH)
    return clip

def padding (src, left=0, top=0, right=0, bottom=0):
    core   = vs.get_core ()
    w      = src.width
    h      = src.height
    clip   = core.fmtc.resample (src, w+left+right, h+top+bottom, -left, -top, w+left+right, h+top+bottom, kernel="point", fulls=True, fulld=True)
    return clip

def cblur16 (src):
    core   = vs.get_core ()
    w      = src.width
    h      = src.height
    blur   = core.fmtc.resample (src, w*4, h*4, kernel="cubic", a1=1, a2=0, fulls=True, fulld=True)
    sharp  = core.fmtc.resample (src, w*4, h*4, kernel="cubic", a1=-1, a2=0, fulls=True, fulld=True)
    dif    = core.std.MakeDiff (blur, sharp)
    dif    = core.fmtc.resample (dif, w, h, kernel="gauss", a1=100, fulls=True, fulld=True)
    clip   = core.std.MergeDiff (src, dif)
    return clip

def max_dif16 (src1, src2, ref):
    core = vs.get_core ()
    clip = core.std.Expr ([src1, src2, ref], ["x z - abs y z - abs > x y ?"])
    return clip

def min_dif16 (src1, src2, ref):
    core = vs.get_core ()
    clip = core.std.Expr ([src1, src2, ref], ["x z - abs y z - abs > y x ?"])
    return clip

def limit_dif16 (flt, src, ref=None, thr=0.25, elast=3.0):
    core  = vs.get_core ()
    thr   = thr * 256
    alpha = 1 / (thr * (elast - 1))
    beta  = elast * thr
    ref   = src if ref is None else ref
    clip  = core.std.Expr ([flt, src, ref], ["x z - abs {thr} <= x x z - abs {beta} >= ? y y {alpha} x y - * {beta} x y - abs - * + ?".format (thr=thr, alpha=alpha, beta=beta)])
    return clip

def merge16_8 (src1, src2, mask):
    core   = vs.get_core ()
    mask16 = core.fmtc.bitdepth (mask, bits=16, fulls=True, fulld=True)
    clip   = core.std.MaskedMerge (src1, src2, mask16)
    return clip

def build_sigmoid_expr (string, inv=False, thr=0.5, cont=6.5):
    core = vs.get_core ()
    x1m0 = "1 {thr} 1 - {cont} * exp 1 + / 1 {cont} {thr} * exp 1 + / -".format (thr=thr, cont=cont)
    x0   = "1 {cont} {thr} * exp 1 + /".format (thr=thr, cont=cont)

    if inv:
       expr = "{thr} 1 " + string + " {x1m0} * {x0} + 0.000001 max / 1 - 0.000001 max log {cont} / -".format (x1m0=x1m0, x0=x0, thr=thr, cont=cont)
    else:
       expr = "1 1 {cont} {thr} " + string + " - * exp + / {x0} - {x1m0} /".format (x1m0=x1m0, x0=x0, thr=thr, cont=cont)
    return expr.format (thr=thr, cont=cont)

def sigmoid_direct (src, thr=0.5, cont=6.5):
    core = vs.get_core ()
    expr = build_sigmoid_expr ("x 65536 /", False, thr, cont)
    clip = core.std.Expr ([src], [expr + " 65536 *"])
    return clip

def sigmoid_inverse (src, thr=0.5, cont=6.5):
    core = vs.get_core ()
    expr = build_sigmoid_expr ("x 65536 /", True, thr, cont)
    clip = core.std.Expr ([src], [expr + " 65536 *"])
    return clip

def clamp16 (src, bright_limit, dark_limit, overshoot=0, undershoot=0):
    core = vs.get_core ()
    os   = overshoot * 256
    us   = undershoot * 256

    clip = core.std.Expr ([src, bright_limit, dark_limit], ["x y {os} + > y {os} + x ? z {us} - < z {us} - x ?".format (os=os, us=us)])
    return clip

def resize16nr (src, w=None, h=None, sx=0, sy=0, sw=0, sh=0, kernel="spline36", kernelh=None, kernelv=None, fh=1, fv=1, taps=4, a1=None, a2=None, a3=None, kovrspl=1, cnorm=True, center=True, fulls=None, fulld=None, cplace="mpeg2", invks=False, invkstaps=4, noring=True):
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
       nrv = round ((1.0 - nrr) * 65535) if nrb else 0
    else:
       nrv = [round ((1.0 - nrr) * 65535), round ((1.0 - nrr) * 65535), round ((1.0 - nrr) * 65535)] if nrb else [0, 0, 0]
    nrm   = core.std.BlankClip (clip=src, width=w, height=h, color=nrv) if nrb and nrf else 0

    main  = core.fmtc.resample (src, w=w, h=h, sx=sx, sy=sy, sw=sw, sh=sh, kernel=kernel, kernelh=kernelh, kernelv=kernelv, fh=fh, fv=fv, taps=taps, a1=a1, a2=a2, a3=a3, kovrspl=kovrspl, cnorm=cnorm, center=center, fulls=fulls, fulld=fulld, cplace=cplace, invks=invks, invkstaps=invkstaps)
    nrng  = core.fmtc.resample (src, w=w, h=h, sx=sx, sy=sy, sw=sw, sh=sh, kernel="gauss", a1=100, center=center, fulls=fulls, fulld=fulld, cplace=cplace) if nrf else main

    clip  = core.rgvs.Repair (main, nrng, 1) if nrf else main
    clip  = core.std.MaskedMerge (main, clip, nrm) if nrf and nrb else clip
    return clip

def ediresample16 (src, w=None, h=None, sx=0, sy=0, sw=None, sh=None, kernel_u="spline64", kernel_d="spline36", taps=4, a1=None, a2=None, a3=None, css=None, fulls=False, fulld=None, cplace="mpeg2", nsize=0, nns=4, qual=2, etype=0, pscrn=1, noring=False, ratiothr=1.0, thr=1.0, elast=1.5):
    core   = vs.get_core ()
    scss   = GetCSS (src)
    css    = scss if css is None else css
    fulld  = fulls if fulld is None else fulld

    Gray   = scss == "GRAY" or css == "GRAY"
    scss   = "444" if Gray else scss
    ocss   = "444" if Gray else css

    InterK = resize16nr if noring else core.fmtc.resample

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
       edgeedi  = EDInter (input, yvct, yhct, 1, 1, True, True, True, nsize, nns, qual, etype, pscrn)
       edgeediU = core.std.ShufflePlanes (edgeedi, planes=1, colorfamily=vs.GRAY)
       edgeediV = core.std.ShufflePlanes (edgeedi, planes=2, colorfamily=vs.GRAY)
       edgeediY = core.std.ShufflePlanes (edgeedi, planes=0, colorfamily=vs.GRAY)
    else:
       edgeediY = core.std.ShufflePlanes (input, planes=0, colorfamily=vs.GRAY) if GetCSS (src) != "GRAY" else input
       edgeediY = edgeediY if Yedit == False else EDInter (edgeediY, yvct, yhct, 1, 1, True, False, False, nsize, nns, qual, etype, pscrn)
       edgeediU = core.std.ShufflePlanes (input, planes=1, colorfamily=vs.GRAY) if Gray == False else 0
       edgeediU = 0 if Uedit == False else EDInter (edgeediU, cvct, chct, 1, 1, Ut, False, False, nsize, nns, qual, etype, pscrn)
       edgeediV = core.std.ShufflePlanes (input, planes=2, colorfamily=vs.GRAY) if Gray == False else 0
       edgeediV = 0 if Vedit == False else EDInter (edgeediV, cvct, chct, 1, 1, Vt, False, False, nsize, nns, qual, etype, pscrn)

    yrh = yrhratio > ratiothr
    yrv = yrvratio > ratiothr
    crh = crhratio > ratiothr
    crv = crvratio > ratiothr

    if enable == False and edit == False:
       0
    else:
       edgeY   = edgeediY if Yedit == False else core.fmtc.resample (edgeediY, ow, oh, (sx * yhrf + yhfix), (sy * yvrf + yvfix), (sw * yhrf), (sh * yvrf), kernelh=(kernel_u if yrh else kernel_d), kernelv=(kernel_u if yrv else kernel_d), taps=taps, a1=a1, a2=a2, a3=a3, fulls=fulls, fulld=fulls)
       edgeU   = 0 if Uedit == False else core.fmtc.resample (edgeediU, owc, ohc, (sxc * chrf + chfix + cphfixe), (syc * cvrf + cvfix), (swc * chrf), (shc * cvrf), kernelh=(kernel_u if crh else kernel_d), kernelv=(kernel_u if crv else kernel_d), taps=taps, a1=a1, a2=a2, a3=a3, fulls=fulls, fulld=fulls)	
       edgeV   = 0 if Vedit == False else core.fmtc.resample (edgeediV, owc, ohc, (sxc * chrf + chfix + cphfixe), (syc * cvrf + cvfix), (swc * chrf), (shc * cvrf), kernelh=(kernel_u if crh else kernel_d), kernelv=(kernel_u if crv else kernel_d), taps=taps, a1=a1, a2=a2, a3=a3, fulls=fulls, fulld=fulls)
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
       flatU = core.std.ShufflePlanes (input, planes=1, colorfamily=vs.GRAY) if (mixed or Uedit == False) and Ut else 0
       flatU = InterK (flatU, owc, ohc, (sxc + cphfix), syc, swc, shc, kernelh=(kernel_u if ch else kernel_d), kernelv=(kernel_u if cv else kernel_d), taps=taps, a1=a1, a2=a2, a3=a3, fulls=fulls, fulld=fulls) if (mixed or Uedit == False) and Ut else 0
       flatV = core.std.ShufflePlanes (input, planes=2, colorfamily=vs.GRAY) if (mixed or Vedit == False) and Vt else 0
       flatV = InterK (flatV, owc, ohc, (sxc + cphfix), syc, swc, shc, kernelh=(kernel_u if ch else kernel_d), kernelv=(kernel_u if cv else kernel_d), taps=taps, a1=a1, a2=a2, a3=a3, fulls=fulls, fulld=fulls) if (mixed or Vedit == False) and Vt else 0	
       flat  = flatY if Gray else core.std.ShufflePlanes ([flatY, flatU, flatV], planes=[0, 0, 0], colorfamily=vs.YUV)

    merge  = (limit_dif16 (flat, edge, thr=thr, elast=elast) if mixed else edge) if edit else flat
    mergeY = core.std.ShufflePlanes (merge, planes=0, colorfamily=vs.GRAY)
    mergeU = core.std.ShufflePlanes (merge, planes=1, colorfamily=vs.GRAY) if Gray == False else 0
    mergeV = core.std.ShufflePlanes (merge, planes=2, colorfamily=vs.GRAY) if Gray == False else 0
    mergeF = core.std.ShufflePlanes ([mergeY, mergeU, mergeV], planes=[0, 0, 0], colorfamily=vs.YUV) if css != "GRAY" else mergeY
    Range  = core.fmtc.bitdepth (mergeF, fulls=fulls, fulld=fulld) if fulls != fulld else mergeF

    return Range

def EDInter (src, vct=1, hct=1, vfield=1, hfield=1, Y=True, U=False, V=False, nsize=0, nns=4, qual=2, etype=0, pscrn=1, honly=False):
    core = vs.get_core ()
    clip = src

    if hct >= 1:
       clip   = clip if honly else core.std.Transpose (clip)
       clip   = core.nnedi3.nnedi3 (clip, hfield, True, Y, U, V, nsize, nns, qual, etype, pscrn)
       hct    = hct - 1
       honly  = hct >= 1
       hfield = 0
       clip   = clip if honly else core.std.Transpose (clip)
    else:
       0

    if vct >= 1 and honly == False:
       clip   = core.nnedi3.nnedi3(clip, vfield, True, Y, U, V, nsize, nns, qual, etype, pscrn)
       vct    = vct - 1
       vfield = 0
    else: 
       0

    clip = (clip if vct <= 0 and hct <= 0 else EDInter (clip, vct, hct, vfield, hfield, Y, U, V, nsize, nns, qual, etype, pscrn, honly)) if Y or U or V else src
    return clip

def GetCSS (src):
    if src.format.id == vs.YUV420P16:
       css = "420"
    elif src.format.id == vs.YUV422P16:
       css = "422"
    elif src.format.id == vs.YUV444P16:
       css = "444"
    elif src.format.id == vs.GRAY16:
       css = "GRAY"
    else:
       css = 0
    return css

def cmsharp16 (src, str=1.0):
    core  = vs.get_core ()
    w     = src.width
    h     = src.height
    super = ediresample16 (src, w*2, h*2, noring=True, kernel_u="spline", kernel_d="spline", taps=3, fulls=True, fulld=True)
    cblur = cblur16 (super)
    cblur = cblur16 (cblur)
    cblur = cblur16 (cblur)
    cdif  = core.std.MakeDiff (cblur, super)
    cdif  = core.fmtc.resample (cdif, w, h, kernel="gauss", a1=100, fulls=True, fulld=True)
    cblur = core.std.MergeDiff (src, cdif)
    mblur = core.ctmf.CTMF (src, radius=1)
    blur  = min_dif16 (cblur, mblur, src)
    dif   = core.std.MakeDiff (src, blur)
    sharp = core.std.MergeDiff (src, dif)
    clip  = sharplimit16 (src, sharp=sharp, str=str)
    return clip

def gensuperclip16 (src, pel=4):
    core   = vs.get_core ()
    w      = src.width
    h      = src.height
    if pel == 2:
       clip = ediresample16 (src, w*2, h*2, sx=0.25, sy=0.25, noring=True, kernel_u="spline", kernel_d="spline", taps=6, fulls=True, fulld=True)
    elif pel == 4:
       clip = ediresample16 (src, w*4, h*4, sx=0.375, sy=0.375, noring=True, kernel_u="spline", kernel_d="spline", taps=6, fulls=True, fulld=True)
    else:
       clip = 0
    return clip

def genlimitclip16 (src, CM=True, str=3.00, repmode=1):
    core   = vs.get_core ()
    if CM:
       sharp = cmsharp16 (src, str=str)
    else:
       sharp = shrinksharp16 (src, str=str)
    rep    = core.rgvs.Repair (sharp, src, repmode)
    clip   = sharplimit16 (src, sharp=rep, str=1.00)
    return clip
