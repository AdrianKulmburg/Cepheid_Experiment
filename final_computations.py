from numpy import *
from uncertainties import ufloat
from uncertainties.umath import *

alpha = ufloat(-2.43, 0.12)
beta = ufloat(4.05, 0.02)

# First for RR Lyr
m_V = ufloat(7.866694890385484, 0.04375793286365991)
P = ufloat(0.3473562451803839, 4.5213693562181337e-07)

M_V = ufloat(0.825, 0.094)

d = 10.0 * pow(10, (m_V - M_V) / 5.0)

print 'For RR Lyr:'
print 'M_V =', M_V
print 'd =', d
print ''

# Then for FF Aql
m_V = ufloat(5.559567492121047, 0.009076871470778491)
P = ufloat(4.552790129452294, 6.777733571442469e-07)

M_V = alpha * (log10(P) - 1.0) - beta

d = 10.0 * pow(10, (m_V - M_V) / 5.0)

print 'For FF Aql:'
print 'M_V =', M_V
print 'd =', d
print ''


# Finally for V 473 Lyr
m_V = ufloat(6.186310206355979, 0.028675321082890122)
P = ufloat(1.5031739836062206, 9.000115360693747e-07)

M_V = alpha * (log10(P) - 1.0) - beta

d = 10.0 * pow(10, (m_V - M_V) / 5.0)

print 'For V 473 Lyr:'
print 'M_V =', M_V
print 'd =', d
print ''

