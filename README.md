# VINEDA -- Volcanic Infrasound Explosions Detector Algorithm

The code implementation of the infrasound detector described in the article. If used, please cite:

``
Vineda: Volcanic Infrasound Explosions Detector Algorithm. Frontiers in Earth Science (2019). Bueno, A., Diaz-Moreno, A., Avarez, I., De la Torre, A., Lamb, O.D., Zuccarello, L. and De Angelis, S.  
``

`Vineda` is available in Python and Matlab.

## Installation & Features (Python)

- Run entirely in Python (2.7). No installation necessary - just use conda or pure Python libraries. 
- Dependencies:
    - Numpy >= 1.10.0
    - Scipy >= 0.19.1
    - Obspy >= 1.1.0

## How to use

A full working example is provided with an infrasound signal. Just call the Python interpreter:

``
python vineda.py
``

Notice that the implementation of the algorithm is defined within ``detect_explosions()`` function, with parameters

- x : the raw infrasound signal (samples).
- fs: original sampling rate (Hz).
- flow: Lower cutoff frequency of the bandpass filter (Hz).
- fhigh: Upper cutoff frequency of the bandpass filter (Hz).
- l: Average duration in seconds of the explosions to be detected.
- nc: Number of sub-bands used in the filter bank.
- beta: Penalty factor for non-impulsive onset

As a result, a tuple containing the CF and the sampling frequency is returned. Vineda can be incorporated
in any workflow by calling this function (with proper params) when needed.  





