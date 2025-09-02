import numpy as np
from matplotlib import pyplot as plt

plot_every = 100

def distance(x1, y1, x2, y2):
    return np.sqrt((x2-x1)**2 + (y2-y1)**2)

def main():

    Nx = 400
    Ny = 100
    tau = .53     # kinematic viscosity/time step
    Nt = 10000

    # lattice speeds and weights

    Nl = 9
    cxs = np.array([0, 0, 1, 1, 1, 0,-1,-1,-1])
    cys = np.array([0, 1, 1, 0,-1,-1,-1, 0, 1])
    weights = np.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36])
    X, Y = np.meshgrid(range(Nx), range(Ny))
    # initial conditions

    F = np.ones([Ny, Nx, Nl]) + .01 * np.random.randn(Ny, Nx, Nl)      # mesoscopic velocity + random perturbations
    F[:, :, 3] = 2.3

    # Cylinder boundary
    cylinder = (X - Nx / 4) ** 2 + (Y - Ny / 2) ** 2 < (Ny / 4) ** 2

    '''
    cylinder = np.full((Ny, Nx), False)    # obstacle shape

    for y in range(0, Ny):
        for x in range(0, Nx):
            if (distance(Nx//4, Ny//2, x, y) < 13 or distance(Nx//4, Ny//2, x, y) < 13):
            #if(distance(0, Ny, x, y)<30 or  distance(0, 0, x, y)<30):
            #or distance(Nx//4, 0, x, y)<30 or distance(Nx//4, Ny, x, y)<30
            #or distance(Nx//2, 0, x, y)<30 or distance(Nx//2, Ny, x, y)<30
            #or distance(Nx*(3/4), 0, x, y)<30 or distance(Nx*(3//4), Ny, x, y)<30):
                cylinder[y][x] = True
    '''


    # main loop

    for it in range(Nt):
        print(it)

        # Zhou - He boundary conditions (absorbtion)

        F[:, -1, [6, 7, 8]] = F[:, -2, [6, 7, 8]]
        F[:, 0, [2, 3, 4]] = F[:, 1, [2, 3, 4]]

        # streaming

        for i, cx, cy in zip(range(Nl), cxs, cys):
            F[:, :, i] = np.roll(F[:, :, i], cx, axis=1)
            F[:, :, i] = np.roll(F[:, :, i], cy, axis=0)

        bndryF = F[cylinder, :]
        bndryF = bndryF[:, [0, 5, 6, 7, 8, 1, 2, 3, 4]]

        # flow properties

        rho = np.sum(F, 2)                   # density (from mass balance)
        ux = np.sum(F * cxs, 2) / rho        # x velocity (from momentum balance)
        uy = np.sum(F * cys, 2) / rho        # y velocity (from momentum balance)

        # boundaries

        F[cylinder, :] = bndryF
        ux[cylinder] = 0
        uy[cylinder] = 0

        # collision

        Feq = np.zeros(F.shape)
        for i, cx, cy, w in zip(range(Nl), cxs, cys, weights):
            Feq[:, :, i] = rho * w * (
                1 + 3 * (cx*ux + cy*uy) + 9 * (cx*ux + cy*uy)**2 / 2 - 3 * (ux**2 + uy**2)/2    # Maxwell - Boltzmann distribution
            )

        F += -(1.0/tau) * (F-Feq)

        # plot

        if (it%plot_every == 0):

            plt.imshow(np.sqrt(ux**2+uy**2))
            plt.pause(0.1)
            plt.cla()

    plt.imshow(np.sqrt(ux**2+uy**2))
    plt.cla()

if __name__=="__main__":
    main()