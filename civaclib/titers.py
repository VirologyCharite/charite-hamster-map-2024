import numpy as np


def logTiterToTiter(logTiter):
    """
    Convert a log titer to titer.

    @param logTiter: The C{float} log titer to convert to a titer.
    """
    return "*" if logTiter == "*" else 2 ** float(logTiter) * 10


def getLogTiter(titer, dilutionStepSize=0):
    """
    Convert a titer to log titer.

    @param titer: the titer to convert
    @param dilutionStepSize: either 1, if discrete titers are used, or 0 if titers
        are continuous.
    """
    if titer == "*":
        return np.nan
    elif titer.startswith("<"):
        return np.log2(int(titer[1:]) / 10) - dilutionStepSize
    elif titer.startswith(">"):
        return np.log2(int(titer[1:]) / 10) + dilutionStepSize
    return np.log2(float(titer) / 10)
