# note: Python syntax here!!!

from functools import reduce


def XOR(l):
    return reduce(lambda x,y: x^y, l)


def AffineSpan(l):
    if len(l) <= 1:
        return l
    off = l[0]
    res = {0}
    for v in l[1:]:
        v ^= off
        if v not in res:
            res |= {u ^ v for u in res}
    res = [u ^ off for u in res]
    return res

def LinearSpan(l):
    if len(l) <= 0:
        return [0]
    res = {0}
    for v in l:
        if v not in res:
            res |= {u ^ v for u in res}
    return res