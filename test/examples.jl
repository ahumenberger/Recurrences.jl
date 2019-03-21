const iv = initvar

cf = Dict()

simpleloop = [
    :(x = 1/2*x)
    :(y = 2y)
]
cf[:simpleloop] = [
    "x(n) = (1/2)^n * $(iv(:x))"
    "y(n) = 2^n * $(iv(:y))"
]

cohencu = [
    :(m = m+1)
    :(x = x+y)
    :(y = y+z)
    :(z = z+6)
]
cf[:cohencu] = []

freire1 = [
    :(x = x-r)
    :(r = r+1)
]
cf[:freire1] = []

freire2 = [
    :(x = x-s)
    :(s = s+6*r+3)
    :(r = r+1)
]
cf[:freire2] = []

euclidex = [[
    :(a = a - b)
    :(p = p - q)
    :(r = r - s)
], [
    :(q = q - p)
    :(b = b - a)
    :(s = s - r)
]]
cf[:(euclidex[1])] = []
cf[:(euclidex[2])] = []

fermat = [[
    :(r = r - v)
    :(v = v + 2)
], [
    :(r = r + u)
    :(u = u + 2)
]]
cf[:(fermat[1])] = []
cf[:(fermat[2])] = []

wensley = [[
    :(b = b/2)
    :(d = d/2)
], [
    :(a = a+b)
    :(y = y+d/2)
    :(b = b/2)
    :(d = d/2)
]]
cf[:(wensley[1])] = []
cf[:(wensley[2])] = []

lcm = [[
    :(x = x - y)
    :(v = v + u)
], [
    :(y = y - x)
    :(u = u + v)
]]
cf[:(lcm[1])] = []
cf[:(lcm[2])] = []

knuth = [[
    :(t  = r)
    :(r  = 2*r-rp+q+d+2)
    :(rp = t)
    :(q  = q+4)
    :(d  = d+2)
], [
    :(t  = r)
    :(r  = 2*r-rp+q)
    :(rp = t)
    :(d  = d+2)
], [
    :(t  = r)
    :(r  = 2*r-rp+q-d-2)
    :(rp = t)
    :(q  = q-4)
    :(d  = d+2)
], [
    :(t  = r)
    :(r  = 2*r-rp+q-2*d-4)
    :(rp = t)
    :(q  = q-8)
    :(d  = d+2)
]]
cf[:(knuth[1])] = []
cf[:(knuth[2])] = []
cf[:(knuth[3])] = []
cf[:(knuth[4])] = []

mannadiv = [[
    :(y1 = y1 + 1)
    :(y2 = 0)
    :(y3 = y3 - 1)
], [
    :(y2 = y2 + 1)
    :(y3 = y3 - 1)
]]
cf[:(mannadiv[1])] = []
cf[:(mannadiv[2])] = []

divbin = [[
    :(x = 2*x)
    :(b = b/2)
], [
    :(b = b/2)
    :(x = 2x+1)
    :(r = r-b)
]]
cf[:(divbin[1])] = []
cf[:(divbin[2])] = []