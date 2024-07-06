import magbottle as mgb

particle_1 = mgb.particle(1, 1, [0], [4], [0], [0], [0.4], [0])
particle_2 = mgb.particle(1, -1, [0], [-4], [0], [0], [0.4], [0])
particle_3 = mgb.particle(1, 1, [0], [0], [0], [0], [0], [0.2])
particle_sist = [particle_1, particle_2, particle_3]
time_sist = mgb.time(0, 50, 0.01)
t = mgb.verlet_velocity(particle_sist, time_sist)

mgb.tray(particle_sist, lim=25)
