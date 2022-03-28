# HRTF Database of HATS from KAIST
Gyeong-Tae Lee, Sang-Min Choi, Byeong-Yun Ko, Yong-Hwa Park

Center for Noise and Vibration Control  
Department of Mechanical Engineering  
Korea Advanced Institute of Science and Technology (KAIST)

## What is HRTF?

The spherical coordinate system and head transverse planes for specifying the location of a sound source are shown in the figure below. The origin of the coordinate system is the center of the head between the entrances to the two ear canals. From the origin, the x, y, and z-axes point to the right ear, front, and top of the head, respectively. The horizontal, median, and lateral planes are defined by these three axes. The position of a sound source is defined in the spherical coordinate system as (r, θ, ϕ). The azimuth θ is the angle between the y-axis and the horizontal projection of the position vector, defined in −180° < θ ≤ +180°, where −90°, 0°, +90°, and +180° indicate the left, front, right, and backward directions respectively in the horizontal plane. The elevation ϕ is the angle between the horizontal plane and the position vector of the sound source, defined in −90° ≤ ϕ ≤ +90°, where −90°, 0°, and +90° represent the bottom, front and top directions respectively in the median plane. The sound emitted from a sound source is diffracted and reflected from the torso, head, and pinna, and then reaches both ears. HRTFs are the acoustic transfer functions due to the sound transmission process that accounts for the overall acoustic filtering effect by human anatomy.

![Fig_01](/images/Fig_01.png)

Fig. 1. Illustrations of (a) a spherical coordinate system and (b) head transverse planes.

## HRTF Measurement System

The HRTF measurement system used in this study is designed to measure HRTFs not only on artificial heads but also on humans. However, since this paper was written to present accurate methods for common issues encountered in the measurement of HRTF database, a standard dummy head, Brüel & Kjær (B&K) HATS Type 4100, is selected as the measurement subject. The HRTF measurement system is installed in the anechoic chamber at Korea Advanced Institute of Science and Technology (KAIST). The size of the KAIST anechoic chamber is 3.6 m in width, 3.6 m in length, and 2.4 m in height, and the cut-off frequency is 100 Hz. Thus, the frequency band of interest of the HRTF measurement system is bounded from 120 Hz to around 20 kHz. In addition, the distance r from the head center of the dummy head to the diaphragm center of a speaker module is set to 1.1 m in consideration of the height of the anechoic chamber. When the distance r exceeds 1 m, the frequency characteristics of HRTFs are not sufficiently affected by the distance change [6], such that the HRTFs become distance-independent and are called far-field HRTFs.  

The range of the sound source azimuth θ is set from −180° to +180°. The range of the sound source elevation ϕ was set from −40° to +90° in consideration of the height of the anechoic chamber. According to several studies [44,45], the angular resolution of a HRTF database should be less than 5° in the horizontal plane and less than 10° in the vertical plane. Therefore, the azimuth and elevation resolutions of the HRTF measurement system are set to 5°, respectively. A turntable driven by a servomotor was designed to realize the azimuth resolution of 5° precisely. As shown in Fig. 3, the measurement subject is mounted on the turntable and rotates so that the speaker array faces the azimuth angle set by the direction-of-arrival (DOA) controller. To realize the 5° elevation resolution, as shown in Fig. 3, the semicircular speaker array composed of speaker modules spaced 5° apart to each other are distributed from −40° to +90°. If an elevation angle is set in the DOA controller, only the switch connected to the corresponding speaker module is turned on in the speaker selector, enabling accurate elevation control of the sound source. The total number of sound source locations where HRTFs are measured is 1,944 (72 points in azimuth × 27 points in elevation).  

As shown in Fig. 3, raw transfer functions, BTFs and OTFs, are measured through the audio interface (Audiomatica CLIO FW-02) connected to the host computer via USB 2.0. The sampling rate of the audio interface is 48 kHz, and the sample size of a measured impulse response is 4,096. In addition, the output signal of the audio interface for full frequency excitation was set as a maximum length sequence (MLS). The MLS from the audio interface is amplified through the speaker amplifier (YBA Heritage A200) and then reproduced as a sound source through the speaker module selected by the speaker selector. The acoustic signal input to the microphone of the dummy head is amplified through the microphone conditioner (B&K NEXUS) and then input to the audio interface. The transfer function is calculated through ensemble average in the audio interface software (Audiomatica CLIO 12 Standard) installed on the host computer, based on the acoustic signals measured eight times.  

![Fig_03](/images/Fig_03.png)

Fig. 3. Block diagram of the HRTF measurement system.



