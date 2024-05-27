import math as m
epsilon = 23.43929111
epsilon = m.radians(epsilon)
k = 0.01720209895
file = open('nebo23.txt', 'r')
ok = file.read().split('\n')
x = float(ok[0]) #прямоугольные координаты
y = float(ok[1])
z = float(ok[2])
Vx = float(ok[3]) #компоненты скорости, отнесенные к плоскости стандартного экватора
Vy = float(ok[4])
Vz = float(ok[5])
t = float(ok[6])
file.close()

def f(x, y): #функция поиска углов из косинусов и синусов
    if y > 0 and x > 0:
        z = m.degrees(m.asin(y))
    elif y > 0 and x < 0:
        z = 180 - m.degrees(m.asin(y))
    elif y < 0 and x < 0:
        z = 180 + abs(m.degrees(m.asin(y)))
    else:
        z = 360 - abs(m.degrees(m.asin(y)))
    return z

# Нахождение начальных данных r,V^2,rr',e,E,M
r = m.sqrt(m.pow(x, 2) + m.pow(y, 2) + m.pow(z, 2)) #расстояние до притягивающего тела (Солнца)
V_2 = m.pow(Vx, 2) + m.pow(Vy, 2) + m.pow(Vz, 2)  #линейная скорость на орбите
a = (m.pow(k, 2) * r) / (2 * m.pow(k, 2) - V_2 * r) #большая полуось орбиты
r_r = x * Vx + y * Vy + z * Vz  #rr'
e = m.sqrt(m.pow((1 - r / a), 2) + m.pow(r_r / (k * m.sqrt(a)), 2))
sinE = r_r / (k * m.sqrt(a) * e)
cosE = (1 - r / a) / e
E = f(cosE, sinE)
E = m.radians(E)
M = E - e * m.sin(E)
n = k/m.pow(a,3/2)

# Нахождение P и Q и их проверка (направляющие косинусы орбиты)
Px = (x / r) * m.cos(E) - (Vx * m.sqrt(a) / k) * m.sin(E)
Py = (y / r) * m.cos(E) - (Vy * m.sqrt(a) / k) * m.sin(E)
Pz = (z / r) * m.cos(E) - (Vz * m.sqrt(a) / k) * m.sin(E)
Qx = x*m.sin(E)/(r*m.sqrt(1-m.pow(e,2)))+Vx*m.sqrt(a)*(m.cos(E)-e)/(k*m.sqrt(1-m.pow(e,2)))
Qy = y*m.sin(E)/(r*m.sqrt(1-m.pow(e,2)))+Vy*m.sqrt(a)*(m.cos(E)-e)/(k*m.sqrt(1-m.pow(e,2)))
Qz = z*m.sin(E)/(r*m.sqrt(1-m.pow(e,2)))+Vz*m.sqrt(a)*(m.cos(E)-e)/(k*m.sqrt(1-m.pow(e,2)))
# Проверка
Sum_P = m.pow(Px,2) + m.pow(Py,2) + m.pow(Pz,2)
Sum_Q=m.pow(Qx,2) + m.pow(Qy,2) + m.pow(Qz,2)
Mult_PQ=Px*Qx+Py*Qy+Pz*Qz
if abs(Sum_Q-1)<0.005 and abs(Sum_P-1)<0.005 and abs(Mult_PQ)<0.005:
    print("Проверка пройдена")
else:
    print("Проверка P и Q не пройдена")
    quit()

# Вычисление i, omega, B_omega
sini=m.sqrt(m.pow((Qz*m.cos(epsilon)-Qy*m.sin(epsilon)),2)+m.pow((Pz*m.cos(epsilon)-Py*m.sin(epsilon)),2))
sin_omega=(Pz*m.cos(epsilon)-Py*m.sin(epsilon))/sini
cos_omega=(Qz*m.cos(epsilon)-Qy*m.sin(epsilon))/sini
omega=f(cos_omega,sin_omega)
sin_Bomega=(Py*m.cos(m.radians(omega))-Qy*m.sin(m.radians(omega)))/m.cos(epsilon)
cos_Bomega=Px*m.cos(m.radians(omega))-Qx*m.sin(m.radians(omega))
Bomega=f(cos_Bomega,sin_Bomega)
cosi=-(Px*m.sin(m.radians(omega))+Qx*m.cos(m.radians(omega)))/m.sin(m.radians(Bomega))
i=f(cosi,sini)
M=m.degrees(M)

# Ответ
rezultat2 = open('itog2.txt', 'w')
rezultat2.write(str(e)+"\n")
rezultat2.write(str(a)+' а.е.'+'\n')
rezultat2.write(str(omega)+'°'+"\n") #аругмент перигелия
rezultat2.write(str(Bomega)+'°'+'\n') #долгота восходящего узла
rezultat2.write(str(i)+'°'+'\n') #наклон пл-ти орбиты
rezultat2.write(str(M)+'°'+'\n') #ср аномалия
rezultat2.close()