# SMS-EPI: Multiband EPI for rodent fMRI
This new multiband-EPI sequence enables simultaneous acquisition of multiple 2D slices regardless of the number of RF receiving channels. It is well-suited for preclinical scanners that don't have high density coil array. It was developed and tested on a Bruker MR scanner (Biospec 94/30) with BGA-12S HP Gradient (max 66G/cm) and a volume Tx coil and a 10mm single-loop Rx coil running ParaVision 6.0.1. It is demonstrated that:

Multiband fMRI

4 slices (4-band) can be simultaneously acquired with increased SNR;
Allow drastically reduction of TR for enhancing temporal resolution;
Enable isotropic resolution by acquiring many thin-slices (eg, 36 slices in 0.5s).

Multiband DTI

4 slices can be simultaneously acquired with slice overlapped;

More than 80 slices with b-value of 3000 s/mm2 can be acquired with TR < 2s;

Enable the use of optimal TR to increase SNR and sensitivity;

Reduce variations of FA, diffusivity and eigenvector;

Reduce CSF contamination.

Key Features:

2 to 4 folds acceleration without array coil.

SliceGRAPPA acceleration with array coil for fMRI.

Both spin-echo and gradient-echo modes are available.

The sequence is now available for academic purposes under a research agreement.

NOTE: no commercial use.

Citation:

For more technical details and citation, please see:

HL Lee, Z Li, EJ Coulson, KH Chuang. “Ultrafast fMRI of the rodent brain using simultaneous multi-slice EPI”, NeuroImage 2019.

HL Lee, X Zhou, Z Li, KH Chuang. “Optimizing Diffusion MRI Acquisition Efficiency of Rodent Brain Using Simultaneous Multislice EPI.” NMR Biomed, 2021. DOI: 10.1002/nbm.4398 

Considerations

The multiband EPI method for fMRI uses longer echo train and extended FOV to avoid slice overlap, which could result in the following issues:

Susceptibility artifact: the sequence is more sensitive to B0 inhomogeneity. Good shimming and/or B0 distortion correction (eg, TOPUP) is recommended.

Resolution: the sequence has larger point spread function in the phase encoding direction. Higher bandwidth and/or faster gradient would be useful.

Multiband RF pulse

Hermite pulse is preferred for low RF power demand but 20% gap is recommended to avoid inter-slice cross-talk.

Sinc7 pulse is for gap-less acquisition but peak power is high, which would need to use longer pulse duration.

Multiband RF pulse needs to match particular slice arrangement. For example, to acquire 5 slices with MB factor = 3 (ie, total = 15 slices) with gap/thickness ratio of 20% needs a pulse like hermite_MB3_s05_r20.exc.

Gradient limit

Very high z-gradient blip will be needed if you want to acquire 2 slices with high MB factor.

The multiband DTI uses Hadamard encoding. To acquire MB=4, the scan needs to repeat 4 times.

White matter T1-optimized TR is recommended. A list of TR at different field strength can be found in the NMR Biomed 2021 paper.

Work-in-progress

multiband ASL

Ultrahigh acceleration (R=8) with array coil

Codes for PV 5.1
