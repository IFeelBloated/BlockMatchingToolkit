import vapoursynth as vs
import math

def ftoi16 (src):
    core   = vs.get_core ()
    clip   = core.fmtc.bitdepth (src, bits=16, flt=False, fulls=True, fulld=True, dmode=1)
    return clip

def i16tof (src):
    core   = vs.get_core ()
    clip   = core.fmtc.bitdepth (src, bits=32, flt=True, fulls=True, fulld=True, dmode=1)
    return clip

def sbrf (src):
    core   = vs.get_core ()
    rg11   = core.rgvs.RemoveGrain (ftoi16 (src), 11)
    rg11   = i16tof (rg11)
    rg11D  = core.std.MakeDiff (src, rg11)
    rg11DR = core.rgvs.RemoveGrain (ftoi16 (rg11D), 11)
    rg11DR = i16tof (rg11DR)
    rg11DD = core.std.Expr ([rg11D, rg11DR], ["x y - x 0.5 - * 0 < 0.5 x y - abs x 0.5 - abs < x y - 0.5 + x ? ?"])
    clip   = core.std.MakeDiff (src, rg11DD)
    return clip

def shrinksharpf (src, str=1.0):
    core   = vs.get_core ()
    blur   = core.std.Median (ftoi16 (src))
    blur   = core.rgvs.Repair (blur, ftoi16 (src), 16)
    blur   = i16tof (blur)
    dif    = core.std.Expr ([src, blur], ["y x - {str} * 3 * 0.5 +".format (str=str)])
    dif    = sbrf (dif)
    difrg  = core.rgvs.RemoveGrain (ftoi16 (dif), 11)
    difrg  = i16tof (difrg)
    fdif   = core.std.Expr ([dif, difrg], ["y 0.5 - abs x 0.5 - abs > y 0.5 ?"])
    clip   = core.std.MergeDiff (src, fdif)
    return clip

def sharplimitf (src, sharp, str=1.0):
    core   = vs.get_core ()
    clip   = core.std.Expr ([src, sharp], ["{x} {y} {x} - abs 4 / log 1 4 / * exp 4 * {str} * {y} {x} - {y} {x} - abs 1.001 + / * + 255 /".format (str=str, x="x 255 *", y="y 255 *")])
    return clip

def gaussf (src, p=30):
    core   = vs.get_core ()
    upsmp  = core.fmtc.resample (src, src.width * 2, src.height * 2, kernel="gauss", a1=100, fulls=True, fulld=True)
    clip   = core.fmtc.resample (upsmp, src.width, src.height, kernel="gauss", a1=p, fulls=True, fulld=True)
    return clip

def xymax (src1, src2):
    core   = vs.get_core ()
    clip   = core.std.Expr ([src1, src2], ["x y > x y ?"])
    return clip

def xymin (src1, src2):
    core   = vs.get_core ()
    clip   = core.std.Expr ([src1, src2], ["x y > y x ?"])
    return clip

def maxmulti (src, start=None, a=2, tr=6):
    core   = vs.get_core ()
    start  = xymax (core.std.SelectEvery (src, tr*2+1, 0), core.std.SelectEvery (src, tr*2+1, 1)) if start is None else start
    max    = xymax (start, core.std.SelectEvery (src, tr*2+1, a))
    a      = a+1
    clip   = max if a == tr*2+1 else maxmulti (src, start=max, tr=tr, a=a)
    return clip

def minmulti (src, start=None, a=2, tr=6):
    core   = vs.get_core ()
    start  = xymin (core.std.SelectEvery (src, tr*2+1, 0), core.std.SelectEvery (src, tr*2+1, 1)) if start is None else start
    min    = xymin (start, core.std.SelectEvery (src, tr*2+1, a))
    a      = a+1
    clip   = min if a == tr*2+1 else minmulti (src, start=min, tr=tr, a=a)
    return clip

def padding (src, left=0, top=0, right=0, bottom=0):
    core   = vs.get_core ()
    w      = src.width
    h      = src.height
    clip   = core.fmtc.resample (src, w+left+right, h+top+bottom, -left, -top, w+left+right, h+top+bottom, kernel="point", fulls=True, fulld=True)
    return clip

def cblurf (src):
    core   = vs.get_core ()
    w      = src.width
    h      = src.height
    blur   = core.fmtc.resample (src, w*4, h*4, kernel="cubic", a1=1, a2=0, fulls=True, fulld=True)
    sharp  = core.fmtc.resample (src, w*4, h*4, kernel="cubic", a1=-1, a2=0, fulls=True, fulld=True)
    dif    = core.std.MakeDiff (blur, sharp)
    dif    = core.fmtc.resample (dif, w, h, kernel="gauss", a1=100, fulls=True, fulld=True)
    clip   = core.std.MergeDiff (src, dif)
    return clip

def max_dif (src1, src2, ref):
    core = vs.get_core ()
    clip = core.std.Expr ([src1, src2, ref], ["x z - abs y z - abs > x y ?"])
    return clip

def min_dif (src1, src2, ref):
    core = vs.get_core ()
    clip = core.std.Expr ([src1, src2, ref], ["x z - abs y z - abs > y x ?"])
    return clip

def limit_dif_f (flt, src, ref=None, thr=0.25, elast=3.0):
    core  = vs.get_core ()
    thr   = thr / 255
    alpha = 1 / (thr * (elast - 1))
    beta  = elast * thr
    ref   = src if ref is None else ref
    clip  = core.std.Expr ([flt, src, ref], ["x z - abs {thr} <= x x z - abs {beta} >= ? y y {alpha} x y - * {beta} x y - abs - * + ?".format (thr=thr, alpha=alpha, beta=beta)])
    return clip

def build_sigmoid_expr (string, inv=False, thr=0.5, cont=6.5):
    x1m0 = "1 {thr} 1 - {cont} * exp 1 + / 1 {cont} {thr} * exp 1 + / -".format (thr=thr, cont=cont)
    x0   = "1 {cont} {thr} * exp 1 + /".format (thr=thr, cont=cont)

    if inv:
       expr = "{thr} 1 " + string + " {x1m0} * {x0} + 0.000001 max / 1 - 0.000001 max log {cont} / -".format (x1m0=x1m0, x0=x0, thr=thr, cont=cont)
    else:
       expr = "1 1 {cont} {thr} " + string + " - * exp + / {x0} - {x1m0} /".format (x1m0=x1m0, x0=x0, thr=thr, cont=cont)
    return expr.format (thr=thr, cont=cont)

def sigmoid_direct_f (src, thr=0.5, cont=6.5):
    core = vs.get_core ()
    expr = build_sigmoid_expr ("x", False, thr, cont)
    clip = core.std.Expr ([src], [expr])
    return clip

def sigmoid_inverse_f (src, thr=0.5, cont=6.5):
    core = vs.get_core ()
    expr = build_sigmoid_expr ("x", True, thr, cont)
    clip = core.std.Expr ([src], [expr])
    return clip

def clampf (src, bright_limit, dark_limit, overshoot=0, undershoot=0):
    core = vs.get_core ()
    os   = overshoot / 255
    us   = undershoot / 255

    clip = core.std.Expr ([src, bright_limit, dark_limit], ["x y {os} + > y {os} + x ? z {us} - < z {us} - x ?".format (os=os, us=us)])
    return clip

def resizefnr (src, w=None, h=None, sx=0, sy=0, sw=0, sh=0, kernel="spline36", kernelh=None, kernelv=None, fh=1, fv=1, taps=4, a1=None, a2=None, a3=None, kovrspl=1, cnorm=True, center=True, fulls=None, fulld=None, cplace="mpeg2", invks=False, invkstaps=4, noring=True):
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

    clip  = i16tof (core.rgvs.Repair (ftoi16 (main), ftoi16 (nrng), 1)) if nrf else main
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

def ediresamplef (src, w=None, h=None, sx=0, sy=0, sw=None, sh=None, kernel_u="spline64", kernel_d="spline36", taps=4, a1=None, a2=None, a3=None, css=None, fulls=False, fulld=None, cplace="mpeg2", nsize=0, nns=4, qual=2, etype=0, pscrn=1, noring=False, ratiothr=1.0, thr=1.0, elast=1.5):
    core   = vs.get_core ()
    scss   = GetCSS (src)
    css    = scss if css is None else css
    fulld  = fulls if fulld is None else fulld

    Gray   = scss == "GRAY" or css == "GRAY"
    scss   = "444" if Gray else scss
    ocss   = "444" if Gray else css

    InterK = resizefnr if noring else core.fmtc.resample

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
    input16 = core.fmtc.bitdepth (input, bits=16, flt=False, fulls=fulls, fulld=fulls)

    if enable == False and edit == False:
       0
    elif yhct == chct and yvct == cvct and scss == "420":
       edgeedi  = EDInter (input16, yvct, yhct, 1, 1, True, True, True, nsize, nns, qual, etype, pscrn)
       edgeediU = core.std.ShufflePlanes (edgeedi, planes=1, colorfamily=vs.GRAY)
       edgeediV = core.std.ShufflePlanes (edgeedi, planes=2, colorfamily=vs.GRAY)
       edgeediY = core.std.ShufflePlanes (edgeedi, planes=0, colorfamily=vs.GRAY)
    else:
       edgeediY = core.std.ShufflePlanes (input16, planes=0, colorfamily=vs.GRAY) if GetCSS (src) != "GRAY" else input16
       edgeediY = edgeediY if Yedit == False else EDInter (edgeediY, yvct, yhct, 1, 1, True, False, False, nsize, nns, qual, etype, pscrn)
       edgeediU = core.std.ShufflePlanes (input16, planes=1, colorfamily=vs.GRAY) if Gray == False else 0
       edgeediU = 0 if Uedit == False else EDInter (edgeediU, cvct, chct, 1, 1, Ut, False, False, nsize, nns, qual, etype, pscrn)
       edgeediV = core.std.ShufflePlanes (input16, planes=2, colorfamily=vs.GRAY) if Gray == False else 0
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
       edge    = core.fmtc.bitdepth (edge, bits=32, flt=True, fulls=fulls, fulld=fulls)

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

    merge  = (limit_dif_f (flat, edge, thr=thr, elast=elast) if mixed else edge) if edit else flat
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

def cmsharpf (src, str=1.0):
    core  = vs.get_core ()
    w     = src.width
    h     = src.height
    super = ediresamplef (src, w*2, h*2, noring=True, kernel_u="spline", kernel_d="spline", taps=3, fulls=True, fulld=True)
    cblur = cblurf (super)
    cblur = cblurf (cblur)
    cblur = cblurf (cblur)
    cdif  = core.std.MakeDiff (cblur, super)
    cdif  = core.fmtc.resample (cdif, w, h, kernel="gauss", a1=100, fulls=True, fulld=True)
    cblur = core.std.MergeDiff (src, cdif)
    mblur = core.std.Median (ftoi16 (src))
    mblur = i16tof (mblur)
    blur  = min_dif (cblur, mblur, src)
    dif   = core.std.MakeDiff (src, blur)
    sharp = core.std.MergeDiff (src, dif)
    clip  = sharplimitf (src, sharp=sharp, str=str)
    return clip

def gensuperclipf (src, pel=4):
    core   = vs.get_core ()
    w      = src.width
    h      = src.height
    if pel == 2:
       clip = ediresamplef (src, w*2, h*2, sx=0.25, sy=0.25, noring=True, kernel_u="spline", kernel_d="spline", taps=6, fulls=True, fulld=True)
    elif pel == 4:
       clip = ediresamplef (src, w*4, h*4, sx=0.375, sy=0.375, noring=True, kernel_u="spline", kernel_d="spline", taps=6, fulls=True, fulld=True)
    else:
       clip = 0
    return ftoi16 (clip)

def genlimitclipf (src, CM=True, str=3.00, repmode=1):
    core   = vs.get_core ()
    if CM:
       sharp = cmsharpf (src, str=str)
    else:
       sharp = shrinksharpf (src, str=str)
    rep    = core.rgvs.Repair (ftoi16 (sharp), ftoi16 (src), repmode)
    rep    = i16tof (rep)
    clip   = sharplimitf (src, sharp=rep, str=1.00)
    return clip

def getvectorsf (src, supersearch, tr=6, pel=4, dct=5, thsad=400):
    core        = vs.get_core ()
    supersoft   = core.mv.Super (ftoi16 (src), pel=pel, chroma=False, hpad=32, vpad=32, pelclip=supersearch, sharp=2, rfilter=4, levels=0)
    supersharp  = core.mv.Super (ftoi16 (src), pel=pel, chroma=False, hpad=32, vpad=32, pelclip=supersearch, sharp=2, rfilter=2, levels=0)
    def search (isb, delta):
        vectors = core.mv.Analyse (supersoft, isb=isb, overlap=16, blksize=32, search=3, chroma=False, truemotion=True, delta=delta, trymany=True, searchparam=16, pelsearch=16, dct=0, levels=0, divide=2, badrange=-24)
        vectors = core.mv.Recalculate (supersoft, vectors, overlap=8, blksize=16, thsad=thsad//2, chroma=False, truemotion=True, search=3, searchparam=16, dct=dct, smooth=1, divide=2)
        vectors = core.mv.Recalculate (supersharp, vectors, overlap=4, blksize=8, thsad=thsad//2, chroma=False, truemotion=True, search=3, searchparam=16, dct=dct, smooth=1, divide=2)
        vectors = core.mv.Recalculate (supersharp, vectors, overlap=2, blksize=4, thsad=thsad//2, chroma=False, truemotion=True, search=3, searchparam=16, dct=dct, smooth=1, divide=0)
        return vectors
    bv          = [search (True, i) for i in range (tr, 0, -1)]
    fv          = [search (False, i) for i in range (1, tr+1)]
    vmulti      = bv + fv
    vmulti      = core.std.Interleave (vmulti)
    return vmulti

def compensatemultif (src, comp, superclip, vmulti, tr=6, pel=4, thsad=400, thscd1=248, thscd2=130):
    core        = vs.get_core ()
    super       = core.mv.Super (ftoi16 (comp), pel=pel, chroma=False, hpad=32, vpad=32, pelclip=superclip, sharp=2, rfilter=2, levels=0)
    def compensate (delta):
        vectors = core.std.SelectEvery (vmulti, tr*2, delta)
        filter  = core.mv.Compensate (ftoi16 (src), super, vectors, thsad=thsad, thscd1=thscd1, thscd2=thscd2)
        return i16tof (filter)
    bcomp       = [compensate (i) for i in range (0, tr)]
    fcomp       = [compensate (i) for i in range (tr, 2*tr)]
    compmulti   = bcomp + [src] + fcomp
    compmulti   = core.std.Interleave (compmulti)
    return compmulti

def degrainnf (src, comp, superclip, vmulti, tr=6, pel=4, thsad=400, thscd1=248, thscd2=130):
    core        = vs.get_core ()
    super       = core.mv.Super (ftoi16 (comp), pel=pel, chroma=False, hpad=32, vpad=32, pelclip=superclip, sharp=2, rfilter=2, levels=0)
    def MDG1 (a):
        bv      = core.std.SelectEvery (vmulti, tr*2, tr-1-a)
        fv      = core.std.SelectEvery (vmulti, tr*2, tr+a)
        MDG     = core.mv.Degrain1 (ftoi16 (src), super, bv, fv, thsad=thsad, thscd1=thscd1, thscd2=thscd2, plane=0, limit=65535)
        return MDG
    MDGMulti    = [MDG1 (i) for i in range (0, tr)]
    MDGMulti    = i16tof (core.std.Interleave (MDGMulti))
    def MDGMerge (start=None, a=2):
        start   = core.std.Merge (core.std.SelectEvery (MDGMulti, tr, 0), core.std.SelectEvery (MDGMulti, tr, 1), 0.5) if start is None else start
        merge   = core.std.Merge (start, core.std.SelectEvery (MDGMulti, tr, a), 1/(a+1))
        a       = a+1
        clip    = merge if a == tr else MDGMerge (start=merge, a=a)
        return clip
    return MDGMerge ()

def writevec (vec, logout):
    core = vs.get_core ()
    w    = vec.get_frame (0).width
    with open (logout, "w") as f:
         print (w, file=f)
    vec  = core.std.CropAbs (vec, width=w, height=1)
    return vec

def readvec (vec, login):
    core = vs.get_core ()
    with open (login, "r") as f:
         w = int (f.readline ())
    vec  = core.raws.Source (vec, w, 1, src_fmt="Y8")
    return vec

def sharpcalmerf (soft, dif, limit, superdif, superlimit, vmulti, login, pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=248, thscd2=130, repmode=13, str=1.00):
    core    = vs.get_core ()
    vmulti  = readvec (vmulti, login)
    blankd  = core.std.Expr ([dif], "0.5")
    MDG     = degrainnf (blankd, dif, superdif, vmulti, tr=tr, pel=pel, thsad=thsadA, thscd1=thscd1, thscd2=thscd2)
    TA      = core.std.MergeDiff (soft, MDG)
    repclip = core.std.MergeDiff (soft, dif)
    TAR     = core.rgvs.Repair (ftoi16 (TA), ftoi16 (repclip), repmode)
    TAR     = i16tof (TAR)
    comp    = compensatemultif (soft, limit, superlimit, vmulti, tr=tr, pel=pel, thsad=thsadL, thscd1=thscd1, thscd2=thscd2)
    max     = maxmulti (comp, tr=tr)
    min     = minmulti (comp, tr=tr)
    TL      = clampf (TAR, max, min, 0, 0)
    SL      = sharplimitf (soft, sharp=TL, str=str)
    return SL

def nlmcleansef (src, local=2.4, nlocal=1.2):
    core    = vs.get_core ()
    lflt    = core.knlm.KNLMeansCL (src, d=0, a=24, s=1, h=local)
    ldif    = core.std.MakeDiff (src, lflt)
    ldif    = ftoi16 (ldif)
    brightd = i16tof (core.std.Inflate (ldif).std.Deflate ())
    darkd   = i16tof (core.std.Deflate (ldif).std.Inflate ())
    ldifp   = core.std.Expr ([brightd, darkd], ["x 0.5 < y x ?"])
    lflt    = core.std.MergeDiff (lflt, ldifp)
    ldif    = core.std.MakeDiff (src, lflt)
    ldifm   = i16tof (core.std.Median (ftoi16 (ldif)))
    ldifp   = core.std.Expr ([ldifm, ldif], ["x 0.5 - y 0.5 - * 0 < 0.5 x 0.5 - abs y 0.5 - abs < x y ? ?"])
    lflt    = core.std.MergeDiff (lflt, ldifp)
    nlflt   = core.knlm.KNLMeansCL (src, d=0, a=24, s=4, h=nlocal)
    clip    = max_dif (nlflt,lflt,src)
    return clip

def spatialnrf (src, NLM, DFT, lowpass=16, safe=1):
    core   = vs.get_core ()
    DFTL   = gaussf (DFT, p=lowpass)
    NLML   = gaussf (NLM, p=lowpass)
    gauss  = gaussf (src, p=safe)
    NLMH   = core.std.MakeDiff (NLM, NLML)
    HNR    = core.std.MergeDiff (DFTL, NLMH)
    HNRL   = gaussf (HNR, p=safe)
    HNRH   = core.std.MakeDiff (HNR, HNRL)
    clip   = core.std.MergeDiff (gauss, HNRH)
    return clip

def difcleansef (spatial, dif, superdif, vmulti, login, pel=4, tr=6, thsad=4800, thscd1=248, thscd2=130, repmode=13):
    core    = vs.get_core ()
    vmulti  = readvec (vmulti, login)
    blankd  = core.std.Expr ([dif], "0.5")
    comp    = degrainnf (blankd, dif, superdif, vmulti, tr=tr, pel=pel, thsad=thsad, thscd1=thscd1, thscd2=thscd2)
    NR      = core.std.MergeDiff (spatial, comp)
    repclip = core.std.MergeDiff (spatial, dif)
    rep     = core.rgvs.Repair (ftoi16 (NR), ftoi16 (repclip), repmode)
    return core.std.MakeDiff (i16tof (rep), spatial)
