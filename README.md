# VINEDA -- Volcanic Infrasound Explosions Detector Algorithm

In this repository, we include the code implementation of a multi-step detection algorithm for infrasonic explosions.
The algorithm is fully described in the following manuscript: 

``
Vineda: Volcanic Infrasound Explosions Detector Algorithm. Frontiers in Earth Science (2019). Bueno, A., Diaz-Moreno, A., Alvarez, I., De la Torre, A., Lamb, O.D., Zuccarello, L. and De Angelis, S.  
``

If used as part of any research work or data processing pipeline, we would appreciate citations. 

## Installation & Features (Python)

- VINEDA runs entirely in Python (2.7 and 3.0) and Matlab 2018. In Python, the library can be directly imported as an 
independent Python package, or interfaced with conda / docker environments.

- Dependencies:
    - Numpy >= 1.10.0
    - Scipy >= 0.19.1
    - Obspy >= 1.1.0

## How to use

A working example is provided with an infrasound signal in ".mat" format. Additionally, from a terminal, just run. 

``
python vineda.py
``

Notice that the implementation of the algorithm is defined within ``detect_explosions()`` function, with the arguments:

- x: A numpy array containing the data samples. From an obspy stream, we can forward a numpy array as stream[0].data
- fs: Original sampling rate of the infrasound signal (Hz)
- flow: Lower cut-off frequency of the bandpass filter (Hz)
- fhigh: Upper cut-off frequency of the bandpass filter (Hz)
- nfb: Number of frequency bands. 
- dmin: Minimum duration of the explosions to be detected (s)
- dmax: Maximum duration of the explosions to be detected (s)
- ndb:  Number of duration bands.
- beta: Factor to reduce non-stationary noises.

As a result, a tuple containing the CF and the sampling frequency is returned. VINEDA can be incorporated
in any data processing workflow by calling this function ``detect_explosions()`` with proper parameter configuration. 
For a practical example, please refer to the jupyter implementation.  
