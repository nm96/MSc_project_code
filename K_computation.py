import simulation

om_1_set = []
B_set = [0.50, 1.65, 1.75]

for B in B_set:
    s = simulation.Simulation()
    s.B = B
    s.T = 2**12
    s.solve()
    s.transform(R=True)
    s.find_peaks()
    om_1_set.append(s.peaks[0])

om_b = om_1_set[0]

K_set = [om_b/om for om in om_1_set]

print(K_set)
