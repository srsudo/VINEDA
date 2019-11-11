import math
import obspy
from obspy.core.trace import Trace, Stats
import numpy as np
import scipy.io
from scipy.signal import butter, lfilter, firwin, filtfilt, medfilt, hilbert
import matplotlib.pyplot as plt


def detect_explosions(x, fs, flow, fhigh, nfb, dmin, dmax, ndb, beta):
    """
    VINEDA main function to detect volcanic explosions from a continuous stream of infrasonic data.

    Args:
        x : (numpy array) A numpy array containing the data samples
        fs: (float) Original sampling rate (Hz)
        flow: (float) Lower cutoff frequency of the bandpass filter (Hz)
        fhigh: (float) Upper cutoff frequency of the bandpass filter (Hz)
        nfb: (int) Number of frequency bands
        dmin: (float) Minimum duration of the explosions to be detected (s)
        dmax: (float) Maximum duration of the explosions to be detected (s)
        ndb: (int) Number of duration bands
        beta: (float) Factor to reduce non-stationary noises

    Returns:
        Numpy Array: Float numpy array with characteristic function
        fsp: float with sampling frequency of the characteristic function
    """

    # Downsampling and decimate the signal

    fsp = math.ceil(2*fhigh)

    # De-trend and decimation
    x.detrend(type='demean')
    data = np.asarray(x.data)
    xr = scipy.signal.resample_poly(data, fsp, fs)

    # Adaptive de-trending filter to remove trends in order to avoid spike noise
    n = int(round(dmax * fsp))

    w_n = float(flow) / float((fsp / 2.0))

    B = firwin(n, w_n)
    trend = filtfilt(B, 1.0, xr)
    xrd = xr - trend

    # Parity by the Kernel of scipy, do check the parity and add 1 otherwise
    med_filt_order = int(round((dmin/10)*fsp))

    if (med_filt_order % 2 ) == 0:
        med_filt_order = (med_filt_order-1)

    #  Median Filter to avoid spiky signals or instrumental noise
    y = medfilt(xrd, med_filt_order)

    # no edge delay window
    fsp = np.round(2*fhigh)

    th = np.arange(0, 1 / float(flow), 1 / fsp)
    haux = np.linspace(1, 0, th.shape[0])

    k = np.arange(1, int(nfb)+1, 1)  # We should keep it one as is freq.

    fl = flow + (k-1)*(fhigh-flow)/float(nfb)
    fh = fl + (fhigh-flow) / float(nfb)

    v = np.zeros((y.shape[0], int(nfb)))
    ei = np.zeros((y.shape[0], int(nfb)))

    for m in range(int(nfb)):
        fc = (fl[m]+fh[m])/2.0
        h1 = haux * np.sin(2 * np.pi * fc * th)
        v[:, m] = lfilter(h1, 1, y)
        ei[:, m] = np.abs(hilbert(v[:, m]))

    e = np.sum(ei, axis=1)

    # Discriminant detector to enhance signals with a sharp rise and gradual decay, whilst mitigating the effect
    # of both, emerging and impulsive noises

    # we create the vector of times
    D = np.linspace(dmin, dmax, ndb)
    di = []

    for dur in D:
        lb = np.round(dur * fsp)

        l1 = np.arange(1, 41, 1)
        l2 = [beta * (x - (lb + 1)) for x in l1]
        h2 = np.concatenate((l1, l2))
        h2 /= np.linalg.norm(h2)

        # CF function building
        cf = lfilter(h2, 1.0, e, axis=0)
        M = int(math.floor(h2.shape[0] / 2.0))  # group delay

        # check this
        cf = np.vstack((cf[M:-1].reshape(-1, 1),
                        np.zeros((M, 1))))

        cf[0:M - 1] = 0
        cf[cf <= 0] = 0
        di.append(cf)

    di = np.hstack(di)
    d = np.prod(di, axis=1)

    # Normalization based on companding method
    cf = companding_mu(d)
    return cf, fsp


def companding_mu(d, mu=1.0, V = float(2**23)):
    """
    Function to apply a non-linear companding method for data normalization whilst emphasizing onset and avoiding lost
    of amplitude information
    Args:
        d: (numpy array) the detected CF of the infrasound stream, without normalization
        mu: (float) the compression factor.
        V: Normalization value. Maximum value of the input signal x (24-bit digitizers)
    Returns:
        cf (numpy array) the normalized cf of infrasound stream.
    """
    cf = 10 * (np.log10(1 + mu * d / V)) / np.log10(1+mu)
    return cf


def read_mat_file(filename):
    """
    Function to read a mat file.
    Args:
        filename: (str) string with the filename we want to read.
    Returns:
        np.array containing the seismic data stream.
    """
    return scipy.io.loadmat(filename)


def plot_stream(infrasound_signal, fsp, graph_title=""):
    """
    Function to plot a infrasonic signal
    Args:
        infrasound_signal: numpy array containing the infrasonic data
        fsp: (float) the sampling frequency of the infrasonic signal
        graph_title: (str) optional parameter with the name of the trace.
    """
    time_vector = np.linspace(0, len(infrasound_signal) / fsp, num=len(infrasound_signal))
    plt.title(graph_title)
    plt.plot(time_vector, infrasound_signal)
    plt.show()
    plt.close()


def convert_to_obspy(signal, fsp):
    # type: (np.array, float) -> obspy.core.trace
    """
    Function to convert a numpy array from any infrasound reading into Obspy format.
    This could help us to take advantage of Obspy optimization routines.
    :param signal: A numpy array containg the signal we want to compute.
    :param fsp: The sampling frequency of the signal.
    :return:
    """
    # we substract the mean

    stats = Stats()
    stats.sampling_rate = float(fsp)
    stats.npts = signal.shape[0]
    return Trace(data=signal.reshape(signal.shape[0], ), header=stats)