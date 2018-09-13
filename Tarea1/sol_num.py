
import math

g = 9.81
NMAX = 10000
TOL = 0.025     #rango de exactitud para los tiempos
skip = 0.05     #incremento para t en checkIntersection

#función F(t) del lado izquierdo de la igualdad en la ecuación trascendente

def F(t):
    return t

'''parte constante del lado derecho G(t) de la igualdad en la ecuación trascendente
el valor theta se debe insertar en grados'''

def constantG(k, vo, theta):
    sinTheta = float(math.sin(math.radians(theta)))
    return float((k * vo * sinTheta + g)/(g * k))

#segunda parte del lado derecho de la iguladad, que al multiplicarse por la constante de constantG, da el lado derecho completo

def G(k, t):
    return (1 - math.exp(-1 * k * t))

'''esta función itera valores de t hasta encontrar una intersección entre las funciones
esta intersección es el tiempo de vuelo
sus argumentos son la constante k y la constante G calculada por constantG'''

def checkIntersection(k, Gcalc):
    t = 1
    while t < NMAX:
        if abs(F(t) - Gcalc * G(k, t)) < TOL:   #si la diferencia entre funciones es menor a la tolerancia, se tiene una respuesta para T
            return float(t)
        else:
            t += skip
        if t == NMAX - 1:
            return 1

#esta función determina el valor de t en el cual se cumple la ecuación trascendental
def solveForT(k, vo, theta):
    Gcalc = float(constantG(k, vo, theta))
    return checkIntersection(k, Gcalc)

#con esta función se calcula la distancia horizontal en un tiempo t
def horizontal(k, vo, theta, t):
    cosTheta = float(math.cos(math.radians(theta)))
    return float((vo * cosTheta)*(1 - math.exp(-1 * k * t))/k)

#con esta función se calcula la distancia vertical en un tiempo t
def vertical(k, vo, theta, t):
    sinTheta = float(math.sin(math.radians(theta)))
    return float((-1 * g * t)/k + (k * vo * sinTheta + g)*(1 - math.exp(-1 * k * t))/(k ** 2))

#con esta función se calcula la distancia horizontal en un tiempo t para k = 0
def horizontalk0(vo, theta, t):
    cosTheta = float(math.cos(math.radians(theta)))
    return float(vo * t * cosTheta)

#con esta función se calcula la distancia vertical en un tiempo t para k = 0
def verticalk0(vo, theta, t):
    sinTheta = float(math.sin(math.radians(theta)))
    return float(((-1 * g * t * t)/2) + (vo * t * sinTheta))

#con esta función se calcula el rango máximo si k = 0
def rangek0(vo, theta):
    sin2Theta = float(math.sin(math.radians(2*theta)))
    return float((vo * vo * sin2Theta)/g)

#con esta función se calcula el tiempo de vuelo si k = 0
def timek0(vo, theta):
    sinTheta = float(math.sin(math.radians(theta)))
    return float((2 * vo * sinTheta)/g)

#con esta función se calcula la velocidad horizontal como la derivada de la función de la posición horizontal
def horizontalVelocity(k, vo, theta, t):
    cosTheta = float(math.cos(math.radians(theta)))
    return float(vo * cosTheta * math.exp(-1 * k * t))

#con esta función se calcula la velocidad vertical como la derivada de la función de la posición vertical
def verticalVelocity(k, vo, theta, t):
    sinTheta = float(math.sin(math.radians(theta)))
    return  float((-1 * g/k) + (((k * vo * sinTheta) + g)* math.exp(-1 * k * t))/k)

#con esta función se calcula la velocidad vertical para el caso especial k = 0
def verticalVelocityk0(vo, theta, t):
    sinTheta = float(math.sin(math.radians(theta)))
    return float((-1 * g * t) + vo * sinTheta)

#con esta función se generan los puntos para graficar distancia horizontal vs vertical
def generate_horizontal_vertical_data(k, vo, theta, flyTime, fileC):
    t = 0

    fileC.write('#k = 0' + '\n')
    fileC.write('x \t y \n')

    while t < timek0(500, 65):      #para el cálculo de k = 0 en el inciso C
        x = horizontalk0(500, 65, t)
        y = verticalk0(500, 65, t)
        fileC.write(str(x) + '\t' + str(y) + '\n')
        t += 1

    t = 0                           #de aquí en adelante es para los cálculos de k != 0

    fileC.write('\n#k = ' + str(k) + '\n')
    fileC.write('x \t y \n')

    while t < flyTime:
        x = horizontal(k, vo, theta, t)
        y = vertical(k, vo, theta, t)
        fileC.write(str(x) + '\t' + str(y) + '\n')
        t += 1

#con esta función se generan los puntos para graficar altura vs tiempo
def generate_height_time_data(k, vo, theta, flyTime, file_height_time):
#se itera t y se calcula la altura en ese tiempo:
    t = 0

    file_height_time.write('#k = 0' + '\n')
    file_height_time.write('t \t y \n')

    while t < timek0(500, 65):      #para los datos de altura vs tiempo con k = 0 en el inciso d
        y = verticalk0(500, 65, t)
        file_height_time.write(str(t) + '\t' + str(y) + '\n')
        t += 1

    t = 0                           #de aquí en adelante es para los cálculos de k != 0

    file_height_time.write('\n#k = ' + str(k) + '\n')
    file_height_time.write('t \t y \n')

    while t < flyTime:
        y = vertical(k, vo, theta, t)
        file_height_time.write(str(t) + '\t' + str(y) + '\n')
        t += 1

#con esta función se generan los puntos para graficar la velocidad horizontal vs tiempo
def generate_horizontal_velocity_time_data(k, vo, theta, flyTime, file_horizontal_velocity_time):
    t = 0

    file_horizontal_velocity_time.write('#k = 0' + '\n')
    file_horizontal_velocity_time.write('t \t vx \n')

    while t < timek0(500, 65):                          #velocidad horizontal para el caso especial k = 0
        vx = horizontalVelocity(0, 500, 65, t)
        file_horizontal_velocity_time.write(str(t) + '\t' + str(vx) + '\n')
        t += 1

    t = 0                                               #velocidad horizontal para k != 0

    file_horizontal_velocity_time.write('\n#k = ' + str(k) + '\n')
    file_horizontal_velocity_time.write('t \t vx \n')

    while t < flyTime:
        vx = horizontalVelocity(k, vo, theta, t)
        file_horizontal_velocity_time.write(str(t) + '\t' + str(vx) + '\n')
        t += 1

#con esta función se generan los puntos para graficar la velocidad horizontal vs tiempo
def generate_vertical_velocity_time_data(k, vo, theta, flyTime, file_vertical_velocity_time):
    t = 0

    file_vertical_velocity_time.write('#k = 0' + '\n')
    file_vertical_velocity_time.write('t \t vy \n')

    while t < timek0(500, 65):                          #velocidad horizontal para el caso especial k = 0
        vy = verticalVelocityk0(500, 65, t)
        file_vertical_velocity_time.write(str(t) + '\t' + str(vy) + '\n')
        t += 1

    t = 0                                               #velocidad horizontal para k != 0

    file_vertical_velocity_time.write('\n#k = ' + str(k) + '\n')
    file_vertical_velocity_time.write('t \t vy \n')

    while t < flyTime:
        vy = verticalVelocity(k, vo, theta, t)
        file_vertical_velocity_time.write(str(t) + '\t' + str(vy) + '\n')
        t += 1

#con esta función se generan los puntos para graficar theta vs rango máximo y se calcula el ángulo que da el mayor alcance
def generate_angle_distance_data(k, vo, file_angle_distance):

    rangeList = []   #aquí se guardarán los valores de los rangos máximos con cada ángulo para k != 0
    rangeListk0 = []    #aquí se guardarán los valores de los rangos máximos con cada ángulo para k = 0
    theta = 0

    file_angle_distance.write('#k = 0\n')
    file_angle_distance.write('theta \t R \n')

    while theta <= 90:
        flyTime = timek0(vo, theta)
        range_i = horizontalk0(500, theta, flyTime)
        rangeListk0.append(range_i)
        file_angle_distance.write(str(theta) + '\t' + str(range_i) + '\n')
        theta += 1

    bestAnglek0 = rangeListk0.index(max(rangeListk0))
    file_angle_distance.write('\nBest angle: ' + str(bestAnglek0) + '\n')

    file_angle_distance.write('\n#k = ' + str(k) + '\n')
    file_angle_distance.write('theta \t R \n')

    theta = 0
    range_i = 0

    while theta <= 90:
        #Gcalc = float(constantG(k, 500, theta))
        #flyTime = checkIntersection(k, Gcalc)
        flyTime = solveForT(k, vo, theta)
        range_i = horizontal(k, 500, theta, flyTime)
        rangeList.append(range_i)
        file_angle_distance.write(str(theta) + '\t' + str(range_i) + '\n')
        theta += 1

    bestAngle = rangeList.index(max(rangeList))
    file_angle_distance.write('\nBest angle: ' + str(bestAngle) + '\n')


#valores en los que se calculará el tiempo de vuelo

#INICIO DE LA EJECUCIÓN

kValues = [0.005, 0.01, 0.05, 0.0025]       #valores usados para el inciso e
#kValues = [0.05, 0.10, 0.15, 0.1]          #valores usados para los demás incisos
voValues = [100, 500, 1000, 1500]
thetaValues = [30, 45, 65, 80]


fileA = open('dataA.txt', 'w')
fileB = open('dataB.txt', 'w')
fileC = open('dataC.txt', 'w')
file_height_time = open('data_height_time.txt', 'w')
file_horizontal_velocity_time = open('horizontal_velocity_time.txt', 'w')
file_vertical_velocity_time = open('vertical_velocity_time.txt', 'w')
file_angle_distance = open('data_angle_distance.txt', 'w')

fileA.write('k \t vo \t theta \t \t G \t \t \t T \n')
fileB.write('k \t R \n')

#en este bloque se calcula el caso especial k = 0 para el inciso B
Tk0, Rk0 = timek0(500, 65), rangek0(500, 65)     #tiempo de vuelo y rango máximo
fileB.write(str(0) + '\t' + str(Rk0) + '\n')


'''en este triple ciclo anidado se calcula T para cada combinación de valores
primero se calcula la constante del lado derecho con constantG
después se revisa la intersección con checkIntersection'''
for k in kValues:
    for vo in voValues:
        for theta in thetaValues:
            flyTime = solveForT(k, vo, theta)
            fileA.write(str(k) + '\t' + str(vo) + '\t' + str(theta) + '\t' + str(flyTime) + '\n')

            if(vo == 500 and theta == 65):           #esto es para el inciso B y el C, donde interesan los valores T con vo = 500 y theta = 65
                maxRange = horizontal(k, vo, theta, flyTime)
                fileB.write(str(k) + '\t' + str(maxRange) + '\n')

                generate_horizontal_vertical_data(k, vo, theta, flyTime, fileC)     #se generan los datos para las gráficas de distancia horizontal vs vertical en dataC
                generate_height_time_data(k, vo, theta, flyTime, file_height_time)
                generate_horizontal_velocity_time_data(k, vo, theta, flyTime, file_horizontal_velocity_time)
                generate_vertical_velocity_time_data(k, vo, theta, flyTime, file_vertical_velocity_time)
                generate_angle_distance_data(k, vo, file_angle_distance)

fileA.close()
fileB.close()
fileC.close()
file_height_time.close()
file_horizontal_velocity_time.close()
file_vertical_velocity_time.close()
file_angle_distance.close()
