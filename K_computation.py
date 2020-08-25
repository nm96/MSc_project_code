import simulation

om_1_set = []
B_set = [1.30, 1.35, 1.40, 1.45, 1.50, 1.55, 1.60, 1.65]

for B in B_set:
    s = simulation.Simulation()
    s.B = B
    s.T = 2**14
    s.solve()
    s.transform(R=True)
    s.find_peaks()
    om_1_set.append(s.peaks[0])

om_b = om_1_set[0]

K_set = [om_b/om for om in om_1_set]

print(K_set)
