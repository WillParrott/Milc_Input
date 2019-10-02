This set us is inteded to produce 2 point, extended source and 3 point MILC input files from a "simplified" yaml file.

The script is currently limited to the case where we have:
1) A 'heavy' parent quark with a variety of masses and spin taste and a single twist.
2) A daughter quark with a single mass and a variety of twists and spin tastes.
3) A spectator quark with one mass and spin and a 'G5-G5' spin taste.

The spectator quark may be the same as one of the daughter quarks, and the spectator quark and any number of the daughter propagators may be loaded from exisiting propagtors.
Additonally, daughter quarks for which we have existing 2 point correlators may be specified. However, please note that if we wish to use these in the three point case, and do not have the propagtors saved, we must regenerate them in the two point case. This does not cost any extra time as we would need to run them fresh in the three point case anyway. 

It is assumed that we will not wish to save parent propagtors. If save is specified, spectator props will be saved, even if they are being loaded. Daughter props will only be saved if save is specified and they are not being loaded. In all cases, if save is not specified, the props will be saved to a temporary file to be read in for the three point calculation and which can be deleted after the three points have run. 
At the present time, only 'G5-G5' daughter propagators will be saved or loaded.

All generated props will be checked in the two point case, and loaded props may be checked or not checked. Whether the propagtor is checked or not will affect the source it is assigned.

No loaded props are checked in the three point case at the current time. 

Calculation set up:

The calculation is set up backwards ie daughter to parent.
This means that the extended source only involves the spin tastes on the parent, and that the cheaper parents are recalculated in the three points, whilst the more expensive daughters are reloaded from temopry files (or from where they were saved if save = True)
