short_f = open(path + 'short.out', 'r')


def coef(g, beta,l_real, l_imag, l_mj, n):
    l = list(range(0,n))
    coef_l = []
    for i in l:
        x = complex(l_real[i], l_imag[i])    # complex
        coef_i = g*beta*(abs(x)**2)*l_mj[i]  # complex
        coef_l.append(coef_i)
    coef = sum(coef_l)
    return coef_l, coef


real_c = []
imag_c = []
mj = []
mi = []
c = 0


for line in short_f:
    vec = line.split()
    if len(vec)<1:
        continue
    elif '#' in line:
        vec2 = vec[2].split(':')
        total_levels = vec2[1]
        continue
    elif '/' in line:
        for i in vec:
            vec2 = i.split('/')
            mj.append(float(vec2[0]))
            mi.append(float(vec2[1]))
    else:
        real_v = vec[:26]
        imag_v = vec[26:]
        for i in real_v:
            real_c.append(float(i))
        for i in imag_v:
            imag_c.append(float(i))


## Constants
#
g = float(7)/6
beta = 0.46686
itot_r = len(real_c)
itot_i = len(imag_c)
i_lvl = itot_r/11
print itot_r, itot_i, i_lvl, 'total levels' + total_levels


coef_lvl1 = coef(g, beta, real_c[0:i_lvl*1], imag_c[0:i_lvl*1], mj[0:i_lvl*1], i_lvl)
coef_lvl2 = coef(g, beta, real_c[(i_lvl*1):(i_lvl*2)], imag_c[(i_lvl*1):(i_lvl*2)], mj[(i_lvl*1):(i_lvl*2)], i_lvl)
coef_lvl3 = coef(g, beta, real_c[(i_lvl*2):(i_lvl*3)], imag_c[(i_lvl*2):(i_lvl*3)], mj[(i_lvl*2):(i_lvl*3)], i_lvl)
coef_lvl4 = coef(g, beta, real_c[(i_lvl*3):(i_lvl*4)], imag_c[(i_lvl*3):(i_lvl*4)], mj[(i_lvl*3):(i_lvl*4)], i_lvl)

coef_lvl5 = coef(g, beta, real_c[(i_lvl*4):(i_lvl*5)], imag_c[(i_lvl*4):(i_lvl*5)], mj[(i_lvl*4):(i_lvl*5)], i_lvl)
coef_lvl6 = coef(g, beta, real_c[(i_lvl*5):(i_lvl*6)], imag_c[(i_lvl*5):(i_lvl*6)], mj[(i_lvl*5):(i_lvl*6)], i_lvl)
coef_lvl7 = coef(g, beta, real_c[(i_lvl*6):(i_lvl*7)], imag_c[(i_lvl*6):(i_lvl*7)], mj[(i_lvl*6):(i_lvl*7)], i_lvl)
coef_lvl8 = coef(g, beta, real_c[(i_lvl*7):(i_lvl*8)], imag_c[(i_lvl*7):(i_lvl*8)], mj[(i_lvl*7):(i_lvl*8)], i_lvl)

coef_lvl9 = coef(g, beta, real_c[(i_lvl*8):(i_lvl*9)], imag_c[(i_lvl*8):(i_lvl*9)], mj[(i_lvl*8):(i_lvl*9)], i_lvl)
coef_lvl10 = coef(g, beta, real_c[(i_lvl*9):(i_lvl*10)], imag_c[(i_lvl*9):(i_lvl*10)], mj[(i_lvl*9):(i_lvl*10)], i_lvl)
coef_lvl11 = coef(g, beta, real_c[(i_lvl*10):(i_lvl*11)], imag_c[(i_lvl*10):(i_lvl*11)], mj[(i_lvl*10):(i_lvl*11)], i_lvl)


d_pend = {}

d_pend['lvl_1'] = coef_lvl1[1]#, coef_lvl1[0]
d_pend['lvl_2'] = coef_lvl2[1]#, coef_lvl1[0]
d_pend['lvl_3'] = coef_lvl3[1]#, coef_lvl1[0]
d_pend['lvl_4'] = coef_lvl4[1]#, coef_lvl1[0]
d_pend['lvl_5'] = coef_lvl5[1]#, coef_lvl1[0]
d_pend['lvl_6'] = coef_lvl6[1]#, coef_lvl1[0]
d_pend['lvl_7'] = coef_lvl7[1]#, coef_lvl1[0]
d_pend['lvl_8'] = coef_lvl8[1]#, coef_lvl1[0]
d_pend['lvl_9'] = coef_lvl9[1]#, coef_lvl1[0]
d_pend['lvl_10'] = coef_lvl10[1]#, coef_lvl1[0]
d_pend['lvl_11'] = coef_lvl11[1]#, coef_lvl1[0]
print d_pend


print coef_lvl1[1]#, coef_lvl1[0]
print coef_lvl2[1]#, coef_lvl2[0]
print coef_lvl3[1]#, coef_lvl3[0]
print coef_lvl4[1]#, coef_lvl4[0]
print coef_lvl5[1]#, coef_lvl4[0]
print coef_lvl6[1]#, coef_lvl4[0]
print coef_lvl7[1]#, coef_lvl4[0]
print coef_lvl8[1]#, coef_lvl4[0]
print coef_lvl9[1]#, coef_lvl4[0]
print coef_lvl10[1]#, coef_lvl1[1]
print coef_lvl11[1]#, coef_lvl1[1]


short_f.close()
