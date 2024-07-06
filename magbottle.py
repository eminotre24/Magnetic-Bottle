import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
from dataclasses import dataclass

#Valores para la simulacion
#Valores de entrada
r_m = 0.1
#Componentes
i_g = np.array([1, 0, 0])
j_g = np.array([0, 1, 0])
k_g = np.array([0, 0, 1])
#Campo Magnetico
k_mf = 10
a = 10

#Dataclasses
@dataclass
class particle:
    mass: float
    charge: float
    xpos: list[float]
    ypos: list[float]
    zpos: list[float]
    xvel: list[float]
    yvel: list[float]
    zvel: list[float]

    #Funciones del objeto
    def get_pos_vector(self, i):
        return np.array([self.xpos[i], self.ypos[i], self.zpos[i]])
    def get_vel_vector(self, i):
        return np.array([self.xvel[i], self.yvel[i], self.zvel[i]])
    def add_pos(self, array):
        self.xpos.append(array[0])
        self.ypos.append(array[1])
        self.zpos.append(array[2])
    def add_vel(self, array):
        self.xvel.append(array[0])
        self.yvel.append(array[1])
        self.zvel.append(array[2])

@dataclass
class time:
    initial: float
    final: float
    deltat: float
      
#Funciones
def gen_matrix_of_arrays(num_arrays, array_size):
    return [np.zeros(array_size) for _ in range(num_arrays)]

Bf = lambda v:(((v[0]+a)*i_g + v[1]*j_g + v[2]*k_g)/((v[0]+a)**2 + v[1]**2 + v[2]**2)**2 - ((v[0]-a)*i_g + v[1]*j_g + v[2]*k_g)/((v[0]-a)**2 + v[1]**2 + v[2]**2)**2)*k_mf*a**(3/2)

def e_f(particles, p, i):
    """Calculo de la aceleracion a una particula "particles[p]" causada por la interaccion electrica
    por el sistema de particulas "particles"

    Args:
        particles (list[particle]): La lista del sistema de particulas. Contiene todas las particulas del sistema.
        p (int): Indice de la particula la cual siente esta aceleracion. La particula esta en la lista particles.
        i (int): Indice temporal para acceder a los atributos de la/s particula/s en un tiempo t(i). No es explicitamente el tiempo.

    Returns:
        numpy.ndarray: Vector de aceleracion en los 3 ejes (x,y,z).
    """
    k = 1
    pp = particles[p].get_pos_vector(i)
    a = np.zeros(len(pp))
    for particle in particles:
        if particle != particles[p]:
            ppe = particle.get_pos_vector(i)
            r = np.linalg.norm(pp - ppe)
            if r < r_m:
                r = r_m
            aa = np.zeros(len(pp))
            for j in range(len(a)):
                aa[j] = (particles[p].charge*particle.charge*k/particles[p].mass)*((pp[j]-ppe[j])/r**3)
            a =+ aa
    return a

def m_ff(particle, vel, i):
    """Calculo de la aceleracion que siente la particula "particle" por el campo magnetico del sistema, es decir,
    el de la botella magnetica.

    Args:
        particle (particle): Particula la cual siente esta aceleracion. Se debe de suministrar el objeto particle.
        vel (numpy.ndarray): Vector en (x, y, z) de velocidad de la particula "particle" en el tiempo en el que se hace esta evaluacion.
        i (int): Indice temporal para acceder a los atributos de la/s particula/s en un tiempo t(i). No es explicitamente el tiempo.

    Returns:
        numpy.ndarray: Vector de aceleracion en los 3 ejes (x,y,z).
    """
    pp = particle.get_pos_vector(i)
    B = Bf(pp)
    a = (particle.charge/particle.mass)*np.cross(vel, B)
    return a

def m_f(particles, p, i, vels=None):
    """Calculo de la aceleracion que siente la particula "particles[p]" por el campo magnetico generado por el movimiento del sistema
    de particulas "particles" (Corrientes).

    Args:
        particles (list[particle]): La lista del sistema de particulas. Contiene todas las particulas del sistema.
        p (int): Indice de la particula la cual siente esta aceleracion. La particula esta en la lista particles.
        i (int): Indice temporal para acceder a los atributos de la/s particula/s en un tiempo t(i). No es explicitamente el tiempo.
        vels (numpy.ndarray, optional): Vector en (x, y, z) de velocidad de la particula "particle" en el tiempo en el que se hace esta evaluacion. Defaults to None.

    Returns:
        numpy.ndarray: Vector de aceleracion en los 3 ejes (x,y,z).
    """
    k = 1 #Constante Biot-Savart
    pp = particles[p].get_pos_vector(i)
    a = np.zeros(len(pp))
    for particle in particles:
        if particle != particles[p]:
            ppe = particle.get_pos_vector(i)
            if vels ==  None:
                vel = particle.get_vel_vector(i)  #Velocidad de la particula que genera el campo
                veli = particles[p].get_vel_vector(i) #Velocidad de la particula que estamos trabajando
            else:
                vel = vels[particles.index(particle)]  #Velocidad de la particula que genera el campo
                veli = vels[p] #Velocidad de la particula que estamos trabajando
            r = np.linalg.norm(pp - ppe)
            if r < r_m:
                r = r_m
            B = (particle.charge*k)*(np.cross(vel,(pp-ppe))/r*3) #Biot-Savart
            aa = (particles[p].charge/particles[p].mass)*np.cross(veli,B)
            a =+ aa
    return a

def verlet_velocity(particles, time):
    """Funcion Velocity Verlet para llevar a cabo la simulacion del sistema de particulas "particles", en un lapso
    de tiempo "time.initial" a "time.final".

    Args:
        particles (list[particle]): La lista del sistema de particulas. Contiene todas las particulas del sistema.
        time (time): Objeto "time" que contiene todos los datos temporales para la simulacion.

    Returns:
        numpy.ndarray: Lista de todos los tiempos en los que se evaluo el sistema.
    """
    h = time.deltat
    t = np.arange(time.initial, time.final, h)
    for i in range(len(t)-1):
        vhalfs = gen_matrix_of_arrays(len(particles), 3)
        for j, particle in enumerate(particles):
            vhalfs[j] = particle.get_vel_vector(i) + 0.5*h*(m_ff(particle, particle.get_vel_vector(i), i) + e_f(particles, j, i) + m_f(particles, j, i))
            newpos = np.zeros(3)
            newpos = particle.get_pos_vector(i) + h*vhalfs[j]
            particles[j].add_pos(newpos)
        for j, particle in enumerate(particles):
            newvel = np.zeros(3)
            newvel = vhalfs[j] + 0.5*h*(m_ff(particle, vhalfs[j], i+1) + e_f(particles, j, i+1) + m_f(particles, j, i+1, vhalfs))
            particles[j].add_vel(newvel)
    return t

def tray(particle_sist, lim):
    """Funcion que despliegua un plot 3D de las trayectorias que siguieron las particulas, junto con el vector field de la botella magnetica, usando la libreria matplotlib.

    Args:
        particles (list[particle]): La lista del sistema de particulas. Contiene todas las particulas del sistema.
        lim (int): Limites del frame de la simulacion
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for particle in particle_sist:
        ax.plot(particle.xpos, particle.ypos, particle.zpos)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_zlabel('Z Position')
    ax.grid(True, color='gray', alpha=0.5, which='both')  # mostramos grid

    #Vector Field del Campo Magnetico
    qv = 20 #Cantidad de vectores
    k = 10
    a = 10

    x = np.linspace(-lim+5, lim-5, qv)
    y = np.linspace(-lim+5, lim-5, qv)
    z = np.linspace(-lim+5, lim-5, qv)
    X, Y, Z = np.meshgrid(x, y, z)

    U = np.zeros_like(X)
    V = np.zeros_like(Y)
    W = np.zeros_like(Z)

    for i in range(len(x)):
        for j in range(len(y)):
            for k in range(len(z)):
                v = np.array([X[i, j, k], Y[i, j, k], Z[i, j, k]])
                B = Bf(v)
                U[i, j, k] = B[0]
                V[i, j, k] = B[1]
                W[i, j, k] = B[2]

    skip = (slice(None, None, 2), slice(None, None, 2), slice(None, None, 2))
    ax.quiver(X[skip], Y[skip], Z[skip], U[skip], V[skip], W[skip], length=3, normalize=True, color='tab:gray', alpha=0.5)

    plt.show()
    
def anim(particle_sist, t, lim):
    """Funcion que despliega una animacion de la trayectoria de las particulas usando la funcion FuncAnimation de la libreria matplotlib.

    Args:
        particles (list[particle]): La lista del sistema de particulas. Contiene todas las particulas del sistema.
        time (time): Objeto "time" que contiene todos los datos temporales para la simulacion.
        lim (int): Limites del frame de la simulacion
    """
    # Definicion de Componentes de la Simulacion
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.grid(True)
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_zlabel('Z Position')

    #Definir la cantidad de objetos de la simulacion de acuerdo a nuestro sistema
    lines = [ax.plot([], [], [], marker='o')[0] for _ in particle_sist]

    #Funcion de definicion del sistema
    def init():
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim, lim)
        ax.grid(True)
        for line in lines:
            line.set_data([], [])
            line.set_3d_properties([])
        return lines

    #Funcion de actualizacion de la animacion
    def update(frame):
        for idx, line in enumerate(lines):
            x_data = particle_sist[idx].xpos[frame]
            y_data = particle_sist[idx].ypos[frame]
            z_data = particle_sist[idx].zpos[frame]
            line.set_data(x_data, y_data)
            line.set_3d_properties(z_data)
        return lines

    ani = FuncAnimation(fig, update, frames=len(t), init_func=init, blit=True, interval=5)
    plt.show()
    