# gallantlab/motion_energy_matlab


motion_evergy_matlab is a software library to calculate motion or static energy features of visual stimuli (movies or images).


## QUICKSTART/USAGE
The codes accept movie/image stimulus inputs and return the energy feature outputs. To test the codes, please prepare the sample stimuli  and run ComputeMotionEnergy.m in MATLAB.

```
>> ComputeMotionEnergy;
```

It allows flexible specifications of filter array parameters, such as the highest spatial frequency, the number of direction of motion, color channels to be used, etc. For example, if you want to load a default parameter set and specify the highest spatial and temporal frequency:

```
>> gparams = preprocWavelets_grid_GetMetaParams(2); % load the default parameter set 2 (as in Nishimoto et al., 2011)
>> gparams.sfmax = 24; % the highest frequency is set to 24 cycles/image
>> gparams.tfmax = 6; % the highest temporal frequency is set to 6 cycles/time window
>> [S_gab, gparams] = preprocWavelets_grid(S_lum, gparams); % calculate the energy features of luminance pattern S_lum
```

For more details, please type ‘help [filename]’ (e.g.,):

>> help preprocWavelets_grid_GetMetaParams


## REQUIREMENTS
-	MATLAB 20xx or higher
-	Image Processing toolbox (for color conversion)


## REFERENCES
The original motion energy model was proposed by Watson and Ahumada (1985) and Adelson and Bergen (1985). The codes were written by Shinji Nishimoto (Nishimoto et al., 2011) and edited and commented by Mark Lescroart.

Watson AB, Ahumada AJ Jr., Model of human visual-motion sensing. J Opt Soc Am A. 2(2):322-41. (1985)
Adelson EH, Bergen JR., Spatiotemporal energy models for the perception of motion., J Opt Soc Am A. 2(2):284-99. (1985)
Nishimoto S, Vu AT, Naselaris T, Benjamini Y, Yu B, Gallant JL., Reconstructing visual experiences from brain activity evoked by natural movies. Curr Biol. 21(19):1641-6. (2011)

