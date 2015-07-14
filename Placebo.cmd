ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Y.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.nlmcleansef (clip, h1=6.4, h2=2.4, h3=1.2) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy NLM.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Y.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.dfttest.DFTTest (clip, sstring="0.0:16.0 0.48:8.0 0.64:0.5 1.0:0.0", sbsize=33, sosize=0, smode=0, tosize=0, tbsize=1, tmode=0) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy DFT.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO Y    = core.raws.Source ("Y.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO NLM  = core.raws.Source ("NLM.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO DFT  = core.raws.Source ("DFT.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.spatialnrf (Y, DFT=DFT, NLM=NLM, lowpass=16, safe=1) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Spatial.rgb -p
del TMP.vpy
del DFT.rgb
del NLM.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Spatial.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.gensuperclipf (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy SuperV.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO clip  = core.raws.Source ("Spatial.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO supv  = core.raws.Source ("SuperV.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO vec   = Placebo.getvectorsf (clip, supv, tr=6, pel=4, dct=5, thsad=400) >> TMP.vpy
ECHO Placebo.writevec (vec, "log.txt").set_output () >> TMP.vpy
call vspipe TMP.vpy Vmulti.rgb -p
del TMP.vpy
del SuperV.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO Y    = core.raws.Source ("Y.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO SPT  = core.raws.Source ("Spatial.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (Y, SPT) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.gensuperclipf (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy SuperD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO spt   = core.raws.Source ("Spatial.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO supd  = core.raws.Source ("SuperD.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO clip  = Placebo.nrfinalf (spt, dif, supd, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsad=4800, thscd1=248, thscd2=130, repmode=13) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy NRFinal.rgb -p
del TMP.vpy
del Dif.rgb
del SuperD.rgb
del Spatial.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("NRFinal.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genlimitclipf (clip, CM=False) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Limit.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Limit.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.gensuperclipf (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy SuperL.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("NRFinal.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (Placebo.shrinksharpf (clip, str=1.00), clip) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.gensuperclipf (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy SuperD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO soft  = core.raws.Source ("NRFinal.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO limit = core.raws.Source ("Limit.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO supd  = core.raws.Source ("SuperD.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO supl  = core.raws.Source ("SuperL.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO clip  = Placebo.sharpcalmerf (soft, dif, limit, supd, supl, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=248, thscd2=130, repmode=13, str=1.00) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy PreS.rgb -p
del TMP.vpy
del Dif.rgb
del Limit.rgb
del SuperD.rgb
del SuperL.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("PreS.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genlimitclipf (clip, CM=True) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Limit.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Limit.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.gensuperclipf (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy SuperL.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("PreS.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (Placebo.cmsharpf (clip, str=1.00), clip) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.gensuperclipf (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy SuperD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO soft  = core.raws.Source ("PreS.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO limit = core.raws.Source ("Limit.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO supd  = core.raws.Source ("SuperD.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO supl  = core.raws.Source ("SuperL.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO clip  = Placebo.sharpcalmerf (soft, dif, limit, supd, supl, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=248, thscd2=130, repmode=13, str=1.00) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy CM.rgb -p
del TMP.vpy
del Dif.rgb
del Limit.rgb
del SuperD.rgb
del SuperL.rgb
del PreS.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("CM.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genlimitclipf (clip, CM=False) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Limit.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Limit.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.gensuperclipf (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy SuperL.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("CM.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (Placebo.shrinksharpf (clip, str=1.00), clip) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.gensuperclipf (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy SuperD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO soft  = core.raws.Source ("CM.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO limit = core.raws.Source ("Limit.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO supd  = core.raws.Source ("SuperD.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO supl  = core.raws.Source ("SuperL.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO clip  = Placebo.sharpcalmerf (soft, dif, limit, supd, supl, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=248, thscd2=130, repmode=13, str=1.00) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy S1.rgb -p
del TMP.vpy
del Dif.rgb
del SuperD.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("S1.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (Placebo.shrinksharpf (clip, str=1.00), clip) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.gensuperclipf (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy SuperD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO soft  = core.raws.Source ("S1.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO limit = core.raws.Source ("Limit.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO supd  = core.raws.Source ("SuperD.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO supl  = core.raws.Source ("SuperL.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO clip  = Placebo.sharpcalmerf (soft, dif, limit, supd, supl, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=248, thscd2=130, repmode=13, str=1.00) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy S2.rgb -p
del TMP.vpy
del Dif.rgb
del SuperD.rgb
del S1.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("S2.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (Placebo.shrinksharpf (clip, str=1.00), clip) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.gensuperclipf (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy SuperD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO soft  = core.raws.Source ("S2.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO limit = core.raws.Source ("Limit.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO supd  = core.raws.Source ("SuperD.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO supl  = core.raws.Source ("SuperL.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO clip  = Placebo.sharpcalmerf (soft, dif, limit, supd, supl, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=248, thscd2=130, repmode=13, str=1.00) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy S3.rgb -p
del TMP.vpy
del Dif.rgb
del SuperD.rgb
del S2.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO CM   = core.raws.Source ("CM.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO S3   = core.raws.Source ("S3.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (S3, CM) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.gensuperclipf (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy SuperD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO soft  = core.raws.Source ("CM.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO limit = core.raws.Source ("Limit.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO supd  = core.raws.Source ("SuperD.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO supl  = core.raws.Source ("SuperL.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO clip  = Placebo.sharpcalmerf (soft, dif, limit, supd, supl, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=248, thscd2=130, repmode=13, str=1.00) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy PreDeconv.rgb -p
del TMP.vpy
del Dif.rgb
del SuperD.rgb
del S3.rgb
del Limit.rgb
del SuperL.rgb
del CM.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("PreDeconv.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.gensuperclipf (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy SuperL.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("NRFinal.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (Placebo.deconvf (clip), clip) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.gensuperclipf (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy SuperD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO soft  = core.raws.Source ("NRFinal.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO limit = core.raws.Source ("PreDeconv.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO supd  = core.raws.Source ("SuperD.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO supl  = core.raws.Source ("SuperL.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO clip  = Placebo.sharpcalmerf (soft, dif, limit, supd, supl, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=248, thscd2=130, repmode=13, str=1.00) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Deconv.rgb -p
del TMP.vpy
del Dif.rgb
del PreDeconv.rgb
del SuperD.rgb
del SuperL.rgb
del NRFinal.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Deconv.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genlimitclipf (clip, CM=False) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Limit.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Limit.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.gensuperclipf (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy SuperL.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Deconv.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (Placebo.shrinksharpf (clip, str=1.00), clip) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.gensuperclipf (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy SuperD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO soft  = core.raws.Source ("Deconv.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO limit = core.raws.Source ("Limit.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO supd  = core.raws.Source ("SuperD.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO supl  = core.raws.Source ("SuperL.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO clip  = Placebo.sharpcalmerf (soft, dif, limit, supd, supl, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=248, thscd2=130, repmode=13, str=1.00) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy S1.rgb -p
del TMP.vpy
del Dif.rgb
del SuperD.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("S1.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (Placebo.shrinksharpf (clip, str=1.00), clip) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.gensuperclipf (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy SuperD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO soft  = core.raws.Source ("S1.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO limit = core.raws.Source ("Limit.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO supd  = core.raws.Source ("SuperD.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO supl  = core.raws.Source ("SuperL.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO clip  = Placebo.sharpcalmerf (soft, dif, limit, supd, supl, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=248, thscd2=130, repmode=13, str=1.00) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy S2.rgb -p
del TMP.vpy
del Dif.rgb
del SuperD.rgb
del S1.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("S2.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (Placebo.shrinksharpf (clip, str=1.00), clip) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.gensuperclipf (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy SuperD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO soft  = core.raws.Source ("S2.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO limit = core.raws.Source ("Limit.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO supd  = core.raws.Source ("SuperD.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO supl  = core.raws.Source ("SuperL.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO clip  = Placebo.sharpcalmerf (soft, dif, limit, supd, supl, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=248, thscd2=130, repmode=13, str=1.00) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy S3.rgb -p
del TMP.vpy
del Dif.rgb
del SuperD.rgb
del S2.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO CM   = core.raws.Source ("Deconv.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO S3   = core.raws.Source ("S3.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (S3, CM) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.gensuperclipf (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy SuperD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO soft  = core.raws.Source ("Deconv.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO limit = core.raws.Source ("Limit.rgb", 736, 480, src_fmt="GRAYS") >> TMP.vpy
ECHO supd  = core.raws.Source ("SuperD.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO supl  = core.raws.Source ("SuperL.rgb", 736 * 4, 480 * 4, src_fmt="GRAY16") >> TMP.vpy
ECHO clip  = Placebo.sharpcalmerf (soft, dif, limit, supd, supl, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=248, thscd2=130, repmode=13, str=1.00) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Retouch.rgb -p
del TMP.vpy
del Dif.rgb
del SuperD.rgb
del S3.rgb
del Limit.rgb
del SuperL.rgb
del Deconv.rgb
