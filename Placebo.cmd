ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Y.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.spatialnr (clip, overshoot=2.4, ring=6.4, nr=3.2, rpass=3, sstring="0.0:16.0 0.48:8.0 0.64:0.5 1.0:0.0", lowpass=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Spatial.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Spatial.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genpelclip (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy PelV.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO clip  = core.raws.Source ("Spatial.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO pelv  = core.raws.Source ("PelV.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO vec   = Placebo.getvectors (clip, pelv, tr=6, pel=4, dct=5, thsad=400) >> TMP.vpy
ECHO Placebo.writevec (vec, "log.txt").set_output () >> TMP.vpy
call vspipe TMP.vpy Vmulti.rgb -p
del TMP.vpy
del PelV.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO Y    = core.raws.Source ("Y.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO SPT  = core.raws.Source ("Spatial.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (Y, SPT) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genpelclip (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy PelD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO spt   = core.raws.Source ("Spatial.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO peld  = core.raws.Source ("PelD.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO clip  = Placebo.nrfinal (spt, dif, peld, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsad=4800, thscd1=10000, thscd2=255, repmode=13) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy NRFinal.rgb -p
del TMP.vpy
del Dif.rgb
del PelD.rgb
del Spatial.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("NRFinal.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genlimitclip (clip, CM=False) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Limit.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Limit.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genpelclip (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy PelL.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("NRFinal.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (Placebo.shrinksharp (clip, str=1.00), clip) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genpelclip (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy PelD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO soft  = core.raws.Source ("NRFinal.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO limit = core.raws.Source ("Limit.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO peld  = core.raws.Source ("PelD.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO pell  = core.raws.Source ("PelL.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO clip  = Placebo.sharpcalmer (soft, dif, limit, peld, pell, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=10000, thscd2=255, repmode=13, str=1.00, safelow=False) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy PreS.rgb -p
del TMP.vpy
del Dif.rgb
del Limit.rgb
del PelD.rgb
del PelL.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("PreS.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genlimitclip (clip, CM=True) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Limit.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Limit.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genpelclip (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy PelL.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("PreS.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (Placebo.cmsharp (clip, str=1.00), clip) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genpelclip (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy PelD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO soft  = core.raws.Source ("PreS.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO limit = core.raws.Source ("Limit.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO peld  = core.raws.Source ("PelD.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO pell  = core.raws.Source ("PelL.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO clip  = Placebo.sharpcalmer (soft, dif, limit, peld, pell, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=10000, thscd2=255, repmode=13, str=1.00, safelow=False) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy CM.rgb -p
del TMP.vpy
del Dif.rgb
del Limit.rgb
del PelD.rgb
del PelL.rgb
del PreS.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("CM.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genlimitclip (clip, CM=False) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Limit.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Limit.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genpelclip (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy PelL.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("CM.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (Placebo.shrinksharp (clip, str=1.00), clip) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genpelclip (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy PelD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO soft  = core.raws.Source ("CM.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO limit = core.raws.Source ("Limit.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO peld  = core.raws.Source ("PelD.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO pell  = core.raws.Source ("PelL.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO clip  = Placebo.sharpcalmer (soft, dif, limit, peld, pell, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=10000, thscd2=255, repmode=13, str=1.00, safelow=False) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy S1.rgb -p
del TMP.vpy
del Dif.rgb
del PelD.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("S1.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (Placebo.shrinksharp (clip, str=1.00), clip) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genpelclip (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy PelD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO soft  = core.raws.Source ("S1.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO limit = core.raws.Source ("Limit.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO peld  = core.raws.Source ("PelD.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO pell  = core.raws.Source ("PelL.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO clip  = Placebo.sharpcalmer (soft, dif, limit, peld, pell, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=10000, thscd2=255, repmode=13, str=1.00, safelow=False) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy S2.rgb -p
del TMP.vpy
del Dif.rgb
del PelD.rgb
del S1.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("S2.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (Placebo.shrinksharp (clip, str=1.00), clip) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genpelclip (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy PelD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO soft  = core.raws.Source ("S2.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO limit = core.raws.Source ("Limit.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO peld  = core.raws.Source ("PelD.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO pell  = core.raws.Source ("PelL.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO clip  = Placebo.sharpcalmer (soft, dif, limit, peld, pell, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=10000, thscd2=255, repmode=13, str=1.00, safelow=False) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy S3.rgb -p
del TMP.vpy
del Dif.rgb
del PelD.rgb
del S2.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO CM   = core.raws.Source ("CM.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO S3   = core.raws.Source ("S3.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (S3, CM) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genpelclip (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy PelD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO soft  = core.raws.Source ("CM.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO limit = core.raws.Source ("Limit.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO peld  = core.raws.Source ("PelD.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO pell  = core.raws.Source ("PelL.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO clip  = Placebo.sharpcalmer (soft, dif, limit, peld, pell, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=10000, thscd2=255, repmode=13, str=1.00, safelow=True) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy PreDeconv.rgb -p
del TMP.vpy
del Dif.rgb
del PelD.rgb
del S3.rgb
del Limit.rgb
del PelL.rgb
del CM.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("PreDeconv.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genpelclip (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy PelL.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("NRFinal.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (Placebo.deconv (clip), clip) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genpelclip (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy PelD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO soft  = core.raws.Source ("NRFinal.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO limit = core.raws.Source ("PreDeconv.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO peld  = core.raws.Source ("PelD.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO pell  = core.raws.Source ("PelL.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO clip  = Placebo.sharpcalmer (soft, dif, limit, peld, pell, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=10000, thscd2=255, repmode=13, str=1.00, safelow=False) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Deconv.rgb -p
del TMP.vpy
del Dif.rgb
del PreDeconv.rgb
del PelD.rgb
del PelL.rgb
del NRFinal.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Deconv.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genlimitclip (clip, CM=False) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Limit.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Limit.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genpelclip (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy PelL.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Deconv.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (Placebo.shrinksharp (clip, str=1.00), clip) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genpelclip (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy PelD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO soft  = core.raws.Source ("Deconv.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO limit = core.raws.Source ("Limit.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO peld  = core.raws.Source ("PelD.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO pell  = core.raws.Source ("PelL.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO clip  = Placebo.sharpcalmer (soft, dif, limit, peld, pell, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=10000, thscd2=255, repmode=13, str=1.00, safelow=False) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy S1.rgb -p
del TMP.vpy
del Dif.rgb
del PelD.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("S1.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (Placebo.shrinksharp (clip, str=1.00), clip) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genpelclip (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy PelD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO soft  = core.raws.Source ("S1.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO limit = core.raws.Source ("Limit.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO peld  = core.raws.Source ("PelD.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO pell  = core.raws.Source ("PelL.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO clip  = Placebo.sharpcalmer (soft, dif, limit, peld, pell, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=10000, thscd2=255, repmode=13, str=1.00, safelow=False) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy S2.rgb -p
del TMP.vpy
del Dif.rgb
del PelD.rgb
del S1.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("S2.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (Placebo.shrinksharp (clip, str=1.00), clip) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genpelclip (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy PelD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO soft  = core.raws.Source ("S2.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO limit = core.raws.Source ("Limit.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO peld  = core.raws.Source ("PelD.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO pell  = core.raws.Source ("PelL.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO clip  = Placebo.sharpcalmer (soft, dif, limit, peld, pell, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=10000, thscd2=255, repmode=13, str=1.00, safelow=False) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy S3.rgb -p
del TMP.vpy
del Dif.rgb
del PelD.rgb
del S2.rgb

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO CM   = core.raws.Source ("Deconv.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO S3   = core.raws.Source ("S3.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = core.std.MakeDiff (S3, CM) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Dif.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core = vs.get_core () >> TMP.vpy
ECHO clip = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO clip = Placebo.genpelclip (clip, pel=4) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy PelD.rgb -p
del TMP.vpy

ECHO import vapoursynth as vs > TMP.vpy
ECHO import Placebo >> TMP.vpy
ECHO core  = vs.get_core () >> TMP.vpy
ECHO soft  = core.raws.Source ("Deconv.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO dif   = core.raws.Source ("Dif.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO limit = core.raws.Source ("Limit.rgb", 704, 576, src_fmt="GRAYS") >> TMP.vpy
ECHO peld  = core.raws.Source ("PelD.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO pell  = core.raws.Source ("PelL.rgb", 704 * 4, 576 * 4, src_fmt="GRAYS") >> TMP.vpy
ECHO clip  = Placebo.sharpcalmer (soft, dif, limit, peld, pell, "Vmulti.rgb", "log.txt", pel=4, tr=6, thsadA=4800, thsadL=400, thscd1=10000, thscd2=255, repmode=13, str=1.00, safelow=True) >> TMP.vpy
ECHO clip.set_output () >> TMP.vpy
call vspipe TMP.vpy Retouch.rgb -p
del TMP.vpy
del Dif.rgb
del PelD.rgb
del S3.rgb
del Limit.rgb
del PelL.rgb
del Deconv.rgb
