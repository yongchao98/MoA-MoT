def hanoi(n, source, target, auxiliary):
    if n == 1:
        print(f"Move disk 1 from Peg {source} to Peg {target}")
    else:
        hanoi(n-1, source, auxiliary, target)
        print(f"Move disk {n} from Peg {source} to Peg {target}")
        hanoi(n-1, auxiliary, target, source)

# Solve the Tower of Hanoi problem for 5 disks
hanoi(5, 1, 2, 3)