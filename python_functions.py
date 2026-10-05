# -*- coding: utf-8 -*-

regF = QuantumRegister(size=nf, name='regF')
regX = QuantumRegister(size=nx, name='regX')
regA = QuantumRegister(size=1, name='regA')

#algorithm Lumen:
def algLumen(circ, xA):
    for i in range (1, nf, 1):
        circ.x(i)
    circ.cx(regA,0)
    listaC=[i for i in range(nq)] #lista qubitow kontrolnych
    for i in range(1, nL, 1):
        j = mojGray(i)
        zmiana = j ^ mojGray((i + 1) % nL)
        omegaDoOmega = np.sqrt(abs(liczKn(nrCor)/liczKn(j)))
        phi_t = 2*phiT*omegaDoOmega
        phi_g = -xA*liczKn(j)-np.pi/2.0
        circ.append(UGate(phi_t, np.pi/2-phi_g, phi_g-np.pi/2).control(nq-1), listaC)
        circ.cx(regA, int(np.log2(zmiana)))#ACNOT
        circ.x(int(np.log2(zmiana)))#ACNOT
    circ.x(regF)

#algorithm U_Int:
def algIntTot(circ):
    for i in range(gray0, ns+gray0, 1):
        if i<ns: ii=i
        else: ii=i-ns
        j = mojGray(ii)
        zmiana = j ^ mojGray((i + 1) % ns)
        x0=(j-ns/2)*dxa #polozenie atomu
        algLumen(circ, x0)
        circ.barrier()
        circ.x(nf+int(np.log2(zmiana)))

