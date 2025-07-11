import sys

# Instruction opcodes
HALT = 0
NOOP = 1
LOAD = 2
STORE = 3
ADD = 4
BITWISENAND = 5
BRANCHIFZERO = 6
BOOP = 7

class BoopMachine:
    """A simulator for the specified toy computer."""
    def __init__(self, R_bits, A_bits, W_bits):
        self.W_bits = W_bits
        self.W_mask = (1 << W_bits) - 1
        
        # We don't need to simulate the full memory/register space,
        # just enough to run the program.
        self.registers = [0] * 256
        self.memory = [0] * 1024
        self.pc = 0
        self.halted = False

    def _to_twos_comp(self, val):
        """Converts a negative number to its two's complement form."""
        if val < 0:
            return ((-val) ^ self.W_mask) + 1
        return val & self.W_mask

    def load_program(self, program, x_input):
        """Initializes memory and registers for a run."""
        for i, instruction in enumerate(program):
            if isinstance(instruction, int):
                self.memory[i] = self._to_twos_comp(instruction)
            else: # It's a tuple representing an instruction
                # This is a simplified encoding for simulation purposes
                self.memory[i] = instruction
        self.registers[0] = x_input

    def run_and_get_stats(self):
        """Runs the loaded program and returns step and boop counts."""
        steps = 0
        boops = 0
        while not self.halted and steps < 100000: # Safety break
            steps += 1
            instruction = self.memory[self.pc]
            op = instruction[0]
            
            # Default behavior is to increment PC
            self.pc += 1

            if op == HALT:
                self.halted = True
            elif op == BOOP:
                boops += 1
            elif op == LOAD: # LOAD reg <- adr
                _, reg_idx, adr = instruction
                self.registers[reg_idx] = self.memory[adr]
            elif op == ADD: # ADD r_dest <- r_src1, r_src2
                _, r_dest, r_src1, r_src2 = instruction
                val = self.registers[r_src1] + self.registers[r_src2]
                self.registers[r_dest] = val & self.W_mask
            elif op == BRANCHIFZERO: # BRANCHIFZERO adr, reg
                _, adr, reg_idx = instruction
                if self.registers[reg_idx] == 0:
                    self.pc = adr
            # Other instructions not needed for this specific program.

        return steps, boops

def analyze_performance():
    """
    Analyzes the performance of the most efficient possible program (a countdown loop)
    and checks if it meets the time complexity requirements.
    """
    
    # Program L: Boops x times using a countdown loop.
    # reg[0] = x (input)
    # reg[1] = -1 (for decrementing)
    # reg[2] = 0 (for unconditional jump)
    program_l = [
        # 0: LOAD reg[1] <- mem[6]   (Load -1 into reg[1])
        (LOAD, 1, 6),
        # 1: BRANCHIFZERO 5, 0       (If x is 0, jump to HALT at addr 5)
        (BRANCHIFZERO, 5, 0),
        # 2: BOOP
        (BOOP,),
        # 3: ADD reg[0] <- reg[0], reg[1] (Decrement x)
        (ADD, 0, 0, 1),
        # 4: BRANCHIFZERO 1, 2       (Unconditional jump to loop start at addr 1, using reg[2] which is always 0)
        (BRANCHIFZERO, 1, 2),
        # 5: HALT
        (HALT,),
        # 6: -1 (data)
        -1
    ]

    print("Analyzing performance of the optimal program (countdown loop):")
    print("-" * 60)
    
    test_values = [0, 1, 10, 15, 16, 17, 20, 50]
    for x in test_values:
        machine = BoopMachine(R_bits=64, A_bits=64, W_bits=512)
        machine.load_program(program_l, x)
        
        steps, boops = machine.run_and_get_stats()
        
        # For x > 0, the steps are 2 (setup) + 4*(x-1) (loops) + 4 (last loop) = 2 + 4x - 4 + 4 = 4x+2.
        # For x = 0, the steps are 2 (LOAD, BRANCH).
        # We can just use the simulator's output.
        
        constraint = x + 50
        passes_test = "PASS" if steps <= constraint else "FAIL"
        
        print(f"Input x = {x:>2} | Boops = {boops:>2} | Steps = {steps:>3} | Constraint (Steps <= {constraint:>3}): {passes_test}")

analyze_performance()