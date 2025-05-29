# Initialize memory and pointer
memory = [0] * 30000
pointer = 0

# Execute the Brainf*ck program
program = "[-]>[-]<>++++++++++[<++++++++++>-]<.++++++++++++++.---..---.<"
pc = 0  # Program counter

while pc < len(program):
    command = program[pc]
    
    if command == '>':
        pointer += 1
    elif command == '<':
        pointer -= 1
    elif command == '+':
        memory[pointer] += 1
    elif command == '-':
        memory[pointer] -= 1
    elif command == '.':
        print(chr(memory[pointer]), end='')
    elif command == '[':
        if memory[pointer] == 0:
            open_brackets = 1
            while open_brackets != 0:
                pc += 1
                if program[pc] == '[':
                    open_brackets += 1
                elif program[pc] == ']':
                    open_brackets -= 1
    elif command == ']':
        if memory[pointer] != 0:
            close_brackets = 1
            while close_brackets != 0:
                pc -= 1
                if program[pc] == '[':
                    close_brackets -= 1
                elif program[pc] == ']':
                    close_brackets += 1
    pc += 1