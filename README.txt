This set up is intended to produce 2 point, extended source and 3 point MILC input files from a "simplified" yaml file.

The script is currently limited to the case where we have:
1) A 'heavy' parent quark with a variety of masses and spin taste and a single twist.
2) A daughter quark with a single mass and a variety of twists and spin tastes.
3) A spectator quark with one mass and spin and a 'G5-G5' spin taste.

The spectator quark may be the same as one of the daughter quarks, and the spectator quark and any number of the daughter propagators may be loaded from exisiting propagtors.
Additonally, daughter quarks for which we have existing 2 point correlators may be specified. However, please note that if we wish to use these in the three point case, and do not have the propagtors saved, we must regenerate them in the two point case. This does not cost any extra time as we would need to run them fresh in the three point case anyway. 

It is assumed that we will not wish to save parent propagtors. If save is specified, spectator props will be saved, even if they are being loaded. Daughter props will only be saved if save is specified and they are not being loaded. In all cases, if save is not specified, the props will be saved to a temporary file to be read in for the three point calculation and which can be deleted after the three points have run. 
At the present time, only 'G5-G5' daughter propagators will be saved or loaded, though other spin tastes will be saved to and loaded from temp. 

All generated props will be checked in the two point case, and loaded props may be checked or not checked. Whether the propagtor is checked or not will affect the source it is assigned.

No loaded props are checked in the three point case at the current time.

My understanding of the check yes option on a loaded propagator is that the propagator is passed to the inverter and the inverter iterates until it reaches the tolerance required (as usual but starting from the loaded propagator). This means that if the propagator you are loading were calculated less precisely than the propagators you are generating, the propagator will be changed slightly. 

Calculation set up:

In the two points:

All loaded check no propagators go on a vector_propagtor_file source.
Loaded check yes propagtors go on a random_color_wall source.
The first generated set of propagators go on a random_color_wall, and all subsequent sets of on the vector_field source reloaded from this.

This is because using a source to generate (and not just check) a propagator means it cannot be reused. However, in order to be saved, a propagtor must be generated from it. This means we must generate the first set of propagators from this source and then reload it for all others. 

In the three points: 

All daughters generated in the two points have been saved to a tempory file so these are loaded in on vector_propagtor_file with check no by default. This means that in the case where the two point for a given daughter propagtor already exists and we are loading in the propagtor for the first time in the three points, it cannot be checked. This is a limitation of the code.  

All parents are recalculated on a vector_field source which is reloaded from the exended source code. 

The calculation is set up backwards (i.e. daughter to parent).
This means that the extended source only involves the spin tastes on the parent, and that the cheaper parents are recalculated in the three points, whilst the more expensive daughters are reloaded from tempory files (or from where they were saved if save = True)


An additional feature is that the tag specified in setting.yaml will be applied to the front of all saved correlators, and also to milc input scripts. The submit script given in settings.yaml is modified to read in tag-milc_2pt_cfg.in and write out tag-milc_2pt_cfg.out etc, and is then resaved as use_this_submit. This is to allow multiple tests on the same configuration to not be saved over each other. Remember to change the name of your submit script in your run-all.sh to "use_this_submit" in order for this to happen. A new submit script will be generated from the original one given each time the code is run. 


To run this code you will need directories and files:

out/
in/setting.yaml
in/input-2pt/
in/input-3pt/
in/input-extsrc/
props/(if you want to save any props)
temp/
submit/ (where submit script will be rewritten to)
