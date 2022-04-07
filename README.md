# HRTF Database of HATS from KAIST
Gyeong-Tae Lee, Sang-Min Choi, Byeong-Yun Ko, Yong-Hwa Park

Center for Noise and Vibration Control  
Department of Mechanical Engineering  
Korea Advanced Institute of Science and Technology (KAIST)

## What is HRTF?

The spherical coordinate system and head transverse planes for specifying the location of a sound source are shown in Fig. 1. In Fig. 1(a), the origin of the coordinate system is the center of the head, between the entrances to the two ear canals. From the origin, the x, y, and z-axes point to the right ear, front, and top of the head, respectively. In Fig. 1(b), the horizontal, median, and lateral planes are defined by these three axes. The position of a sound source is defined in the spherical coordinate system as (r, θ, ϕ). The azimuth θ is the angle between the y-axis and the horizontal projection of the position vector, defined as −180° < θ ≤ +180°, where −90°, 0°, +90°, and +180° indicate the left, front, right, and backward directions, respectively, on the horizontal plane. The elevation ϕ is the angle between the horizontal plane and the position vector of the sound source, defined as −90° ≤ ϕ ≤ +90°, where −90°, 0°, and +90° represent the bottom, front, and top directions, respectively, in the median plane.

![Fig01](/images/Fig01.png)

Fig. 1. Illustrations of (a) spherical coordinate system and (b) head transverse planes.

The sound emitted from a sound source is diffracted and reflected from the torso, head, and pinna, and then reaches both ears as shown in Fig. 2. 

![Fig02](/images/Fig02.png)

Fig. 2. Diagrams of sound transmission from sound source to both ears: (a) top view; (b) side view.

Head-related transfer functions (HRTFs) are acoustic transfer functions due to the sound transmission process that account for the overall acoustic filtering effect by human anatomy. A far-field HRTF of the left or right ear for a sound source of P<sub>S</sub> (r,θ,ϕ) is defined as follows:

![Eq01](/equations/Eq01.png)

where P<sub>L,R</sub> is the complex-valued sound pressure in the frequency domain at the entrance of the left or right ear canal of a subject; P<sub>0</sub> is the complex-valued sound pressure in the frequency domain at the center of the subject’s head in the absence of the subject; the subscripts L and R denote the left and right ears, f refers to frequency, and s refers to a set of parameters related to the dimensions of the subject’s anatomical structures. Although P<sub>L,R</sub> and P<sub>0</sub> are functions of distance r, the effects of r on P<sub>L,R</sub> and P<sub>0</sub> can be regarded as identical under the far-field assumption, so that the effects of r can be canceled out in H<sub>L,R</sub>. Even though Eq. (1) is expressed in terms of ideal sound pressures, when actually measuring HRTFs, the transfer function between the measured sound pressure and the input signal from the measurement system is used. Therefore, it is useful to express HRTFs based on the measured transfer functions. Regarding this point, Eq. (1) can be re-written as follows:

![Eq02](/equations/Eq02.png)

where P ̃<sub>L,R</sub> is the measured sound pressure at the entrance of the left or right ear, P ̃<sub>0</sub> is the measured sound pressure at the head center, and X is the input signal. Both measurement values P ̃<sub>L,R</sub> and P ̃<sub>0</sub> are determined by P<sub>L,R</sub> and P<sub>0</sub>, as well as by H<sub>s</sub> (ϕ,f), which is the transfer function of the measurement system consisting of digital-to-analog converter (DAC), speaker amplifier, speaker module at elevation ϕ in a vertical speaker array, microphone, microphone conditioner, and analog-to-digital converter (ADC). The numerator of Eq. (2) is the binaural transfer function (BTF), which denotes the transfer function between the measured sound pressure at the left or right ear and the input signal; the denominator is the origin transfer function (OTF), which denotes the transfer function between the measured sound pressure at the head center and the input signal. The BTF and OTF are respectively defined as follows:

![Eq03-04](/equations/Eq03-04.png)

In general, when measuring far-field HRTFs, both s and r are constant because the measurement subject is predetermined and the distance of a speaker module from the head center is also fixed for a specific measurement setup. Therefore, based on the BTF and OTF to be measured, the far-field HRTF is defined as follows:

![Eq05](/equations/Eq05.png)

## HRTF Measurement System

The HRTF measurement system used in this study is designed to measure HRTFs not only on artificial heads but also on humans. However, since this work was done to present accurate methods for common issues encountered in the measurement of HRTFs, a standard dummy head, Brüel & Kjær (B&K) HATS Type 4100, is selected as the measurement subject. The HRTF measurement system is installed in the anechoic chamber at the Korea Advanced Institute of Science and Technology (KAIST). The size of the KAIST anechoic chamber is 3.6 m in width, 3.6 m in length, and 2.4 m in height, and the cut-off frequency is 100 Hz. Thus, the frequency band of interest of the HRTF measurement system is bounded from 120 Hz to around 20 kHz. In addition, the distance r from the head center of the dummy head to the diaphragm center of the speaker module is set to 1.1 m in consideration of the height of the anechoic chamber. When the distance r exceeds 1 m, the frequency characteristics of HRTFs are not sufficiently affected by the distance change, such that the HRTFs become distance-independent and are called far-field HRTFs.

The range of the sound source azimuth θ is set from −180° to +180°. The range of the sound source elevation ϕ was set from −40° to +90° in consideration of the height of the anechoic chamber. According to several studies, the angular resolution of an HRTF database should be less than 5° in the horizontal plane and less than 10° in the vertical plane. Therefore, both the azimuth and elevation resolutions of the HRTF measurement system were set to 5°. A turntable driven by a servomotor was designed to precisely realize the azimuth resolution of 5°. As shown in Fig. 3, the measurement subject was mounted on the turntable and rotated so that the speaker array faced the azimuth angle set by the direction-of-arrival (DOA) controller. To realize the 5° elevation resolution, as shown in Fig. 3, a semicircular speaker array was formed using speaker modules spaced 5° apart from each other and distributed from −40° to +90°. Elevation angle was set in the DOA controller, so only the switch connected to the corresponding speaker module was turned on in the speaker selector, enabling accurate elevation control of the sound source. The total number of sound source locations where HRTFs were measured was 1,944 (72 points in azimuth × 27 points in elevation).

As shown in Fig. 3, the raw transfer functions, the BTFs and OTFs, were measured through the audio interface (Audiomatica CLIO FW-02) connected to the host computer via USB 2.0. The sampling rate of the audio interface was 48 kHz; the sample size of the measured impulse responses was 4,096. In addition, the output signal of the audio interface for full frequency excitation was set as a maximum length sequence (MLS). The MLS from the audio interface was amplified through the speaker amplifier (YBA Heritage A200) and then reproduced as a sound source through the speaker module selected by the speaker selector. The acoustic signal input to the microphone of the dummy head was amplified through the microphone conditioner (B&K NEXUS) and then input to the audio interface. The transfer function was calculated using the ensemble average in the audio interface software (Audiomatica CLIO 12 Standard) installed on the host computer, based on the acoustic signals measured eight times.

![Fig03](/images/Fig03.png)

Fig. 3. Block diagram of HRTF measurement system.

![Fig04](/images/Fig04.png)

Fig. 4. HRTF measurement setup with speaker array and dummy head on rotating turntable in anechoic chamber at KAIST: (a) dummy head positioning; (b) speaker array positioning; (c) installation completed in the anechoic chamber; and (d) floor finished with sound-absorbing material.

## HRTF Measurement according to Azimuth and Elevation

The time and frequency domain characteristics of the derived HRTFs were presented as contour maps for each major azimuth and elevation angle. Figs. 5–7 show left and right HRTF pairs as impulse responses, magnitude responses in dB scale, and phase responses in radians for the entire elevation range when the azimuth angles are 0°, 45°, and 90°, respectively.

![Fig05](/images/Fig05.png)

Fig. 5. HRTFs at θ=0˚, H<sub>L,R</sub> (0,ϕ,f): (a1) left HRIRs; (a2) right HRIRs; (b1) magnitude of left HRTFs; (b2) magnitude of right HRTFs; (c1) phase of left HRTFs; and (c2) phase of right HRTFs.

![Fig06](/images/Fig06.png)

Fig. 6. HRTFs at θ=45˚, H<sub>L,R</sub> (45,ϕ,f): (a1) left HRIRs; (a2) right HRIRs; (b1) magnitude of left HRTFs; (b2) magnitude of right HRTFs; (c1) phase of left HRTFs; and (c2) phase of right HRTFs.

![Fig07](/images/Fig07.png)

Fig. 7. HRTFs at θ=90˚, H<sub>L,R</sub> (90,ϕ,f): (a1) left HRIRs; (a2) right HRIRs; (b1) magnitude of left HRTFs; (b2) magnitude of right HRTFs; (c1) phase of left HRTFs; and (c2) phase of right HRTFs.

Figs. 8–10 show left and right HRTF pairs as impulse responses, magnitude responses in dB scale, and phase responses in radians for the entire azimuth range when the elevation angles are 0°, 60°, and 90°, respectively. 

![Fig08](/images/Fig08.png)

Fig. 8. HRTFs at ϕ=0˚, H<sub>L,R</sub> (θ,0,f): (a1) left HRIRs; (a2) right HRIRs; (b1) magnitude of left HRTFs; (b2) magnitude of right HRTFs; (c1) phase of left HRTFs; and (c2) phase of right HRTS.

![Fig09](/images/Fig09.png)

Fig. 9. HRTFs at ϕ=60˚, H<sub>L,R</sub> (θ,60,f): (a1) left HRIRs; (a2) right HRIRs; (b1) magnitude of left HRTFs; (b2) magnitude of right HRTFs; (c1) phase of left HRTFs; and (c2) phase of right HRTFs.

![Fig10](/images/Fig10.png)

Fig. 10. HRTFs at ϕ=90˚, H<sub>L,R</sub> (θ,90,f): (a1) left HRIRs; (a2) right HRIRs; (b1) magnitude of left HRTFs; (b2) magnitude of right HRTFs; (c1) phase of left HRTFs; and (c2) phase of right HRTFs.

## Description of Folders

The files to download are BTF and OTF measurement results, source codes for building HRTF database, and data files about derived HRTFs and binaural sound localization cues, etc.

* The `00_Data` folder consists of HRTF database (`HRTF_HATS`), ILD data (`ILD_HATS`), wideband ILD data (`ILD_HATS_full`), ITD data (`ITD_HATS`), spectral cue data (`PN_HATS`), BTF measurements (`TF_HATS`), and OTF measurements (`TF_Ref`).

* The `01_TF_check` folder consists of scripts to visualize the BTF and OTF measurement results.

* The `02_HRTF_gen` folder consists of scripts to generate the HRTF database.

* The `03_HRTF_check` folder consists of scripts to visualize the HRTF database.

* The `04_ITD-ILD` folder consists of scripts to generate and visualize the ITD, ILD, and wideband ILD data.

* The `05_Directivity` folder consists of scripts to visualize horizontal plane directivity (HPD).

* The `06_PN` folder consists of scripts to generate and visualize the spectral cue (SC) data.

## Citation

If you want to know more about the system design, measurement, and post-processing of the HRTF database, please refer to the paper below.

Also, if you use the presented HRTF database, please cite the paper below.

[Gyeong-Tae Lee, Sang-Min Choi, Byeong-Yun Ko, and Yong-Hwa Park, “HRTF measurement for accurate sound localization cues”, ArXiv Preprint, ArXiv:2203.03166, 2022.](http://arxiv.org/abs/2203.03166)

## License

Except for the `dirplot_F.m` script in the `05_Directivity` folder, the rest of the repository is licensed under the [KAIST License](LICENSE.md).
