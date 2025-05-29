from ortools.sat.python import cp_model

def main():
    # Create the model
    model = cp_model.CpModel()

    # Define the numbers
    numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

    # Create variables for each letter
    A = model.NewIntVarFromDomain(cp_model.Domain.FromValues(numbers), 'A')
    B = model.NewIntVarFromDomain(cp_model.Domain.FromValues(numbers), 'B')
    C = model.NewIntVarFromDomain(cp_model.Domain.FromValues(numbers), 'C')
    D = model.NewIntVarFromDomain(cp_model.Domain.FromValues(numbers), 'D')
    E = model.NewIntVarFromDomain(cp_model.Domain.FromValues(numbers), 'E')
    F = model.NewIntVarFromDomain(cp_model.Domain.FromValues(numbers), 'F')
    G = model.NewIntVarFromDomain(cp_model.Domain.FromValues(numbers), 'G')
    H = model.NewIntVarFromDomain(cp_model.Domain.FromValues(numbers), 'H')
    I = model.NewIntVarFromDomain(cp_model.Domain.FromValues(numbers), 'I')
    J = model.NewIntVarFromDomain(cp_model.Domain.FromValues(numbers), 'J')
    K = model.NewIntVarFromDomain(cp_model.Domain.FromValues(numbers), 'K')
    L = model.NewIntVarFromDomain(cp_model.Domain.FromValues(numbers), 'L')
    M = model.NewIntVarFromDomain(cp_model.Domain.FromValues(numbers), 'M')

    # Add constraints
    model.Add(E > C)
    model.Add(G == 4.8 * M)
    model.Add(L - C == 12)
    model.Add(D + H == 66)
    model.Add(L - D == -35)
    model.Add(E == 3.6 * B)
    model.Add(B - D == -40)
    model.Add(H - G == -8)
    model.Add(C == 1.5 * K)
    model.Add(H == 1.6 * B)
    model.Add(J - E == 9)

    # All variables must take different values
    model.AddAllDifferent([A, B, C, D, E, F, G, H, I, J, K, L, M])

    # Create the solver and solve
    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    if status == cp_model.FEASIBLE or status == cp_model.OPTIMAL:
        result = [solver.Value(A), solver.Value(B), solver.Value(C), solver.Value(D),
                  solver.Value(E), solver.Value(F), solver.Value(G), solver.Value(H),
                  solver.Value(I), solver.Value(J), solver.Value(K), solver.Value(L),
                  solver.Value(M)]
        print(f"<<<{result}>>>")
    else:
        print("No solution found.")

if __name__ == '__main__':
    main()