import math

def raspy_sum(a, b):
    return [x + y for x, y in zip(a, b)]

def raspy_mod(a, val):
    return [x % val for x in a]

class RaspyRunner:
    def __init__(self, tokens):
        self.tokens = list(tokens)
        self.length = len(self.tokens)
        self.indices = list(range(self.length))

    def where(self, cond, true_val, false_val):
        if not isinstance(true_val, list):
            true_val = [true_val] * self.length
        if not isinstance(false_val, list):
            false_val = [false_val] * self.length
        return [t if c else f for c, t, f in zip(cond, true_val, false_val)]

    def a(self, seq):
        return [ord(x) - ord('0') if x.isdigit() else 0 for x in seq]

    def c(self, seq):
        return "".join([chr(x) for x in seq])

    def d(self, seq):
        res = []
        current_sum = 0
        for x in seq:
            current_sum += (1 if x else 0)
            res.append(current_sum)
        return res

    def f(self, i, default, seq):
        return [seq[q - i] if 0 <= q - i < self.length else default for q in self.indices]

    def h(self, i, default, seq):
        return [seq[q + i - 1] if 0 <= q + i - 1 < self.length else default for q in self.indices]
    
    def i_func(self, i, default, seq):
        x = [seq[q - i + 3] if 0 <= q - i + 3 < self.length else default for q in self.indices]
        x = [x[q + i - 3] if 0 <= q + i - 3 < self.length else default for q in self.indices]
        return x

    def j(self, seq):
        min_val = 9999
        for x in seq:
            if isinstance(x, int) or isinstance(x, float):
                min_val = min(min_val, x)
        return [min_val] * self.length

    def l(self, default, sop):
        c = sop.count("_")
        return [sop[i - c] if 0 <= i - c < self.length else default for i in self.indices]

    def m(self, v, i, default="0"):
        sop_str = "".join(self.tokens)
        try:
            split_point = sop_str.index(v)
        except ValueError:
            split_point = -1
        
        if i:
            seq_i_list = [self.tokens[idx] if idx < split_point else "_" for idx in self.indices]
            seq_i = "".join(seq_i_list)
            return self.l(default, seq_i)
        else:
            return [self.tokens[idx] if idx > split_point else default for idx in self.indices]

    def n(self, match, seq):
        y = list(seq)
        for i in range(self.length):
            if not match[i]:
                # Find next j where match is true
                next_j_val = None
                for j in range(i + 1, self.length):
                    if match[j]:
                        next_j_val = seq[j]
                        break
                if next_j_val is not None:
                    y[i] = next_j_val
        return self.where(match, seq, y)

    def s(self, sop):
        # Simplified based on analysis: returns 1 if '7' is in sop, else 0.
        has_7 = False
        for char in sop:
            if char == '7':
                has_7 = True
                break
        
        a_counter = self.where([c == '7' for c in sop], 1, 0)
        a_sum = self.d(a_counter)
        last_index = self.where([i > 1 for i in self.indices], a_sum, "_")
        
        all_last_index_val = 0
        if self.length > 2 and last_index[2] != "_":
            all_last_index_val = last_index[2]

        return [all_last_index_val] * self.length
        
    def q(self, default="_"):
        return self.where([i < 3 for i in self.indices], self.tokens, default)
        
    def r(self, default="_"):
        return self.where([(i > 2) and (i < 6) for i in self.indices], self.tokens, default)

    def p(self, default="_"):
        return self.where([i > self.length - 4 for i in self.indices], self.tokens, default)

    def t(self, seq):
        first_good_idx = self.length
        for i, char in enumerate(seq):
            if char != "_":
                first_good_idx = i
                break
        
        shifted = [seq[q + first_good_idx] if 0 <= q + first_good_idx < self.length else "_" for q in self.indices]
        return shifted


    def u(self):
        # Simplified based on analysis: checks for '7' in three chunks
        aa_chunk = self.q()
        bb_chunk = self.t(self.r())
        cc_chunk = self.t(self.p())

        s_aa = self.s(aa_chunk)
        s_bb = self.s(bb_chunk)
        s_cc = self.s(cc_chunk)

        h_s_aa = self.h(self.length, 0, s_aa)
        h_s_bb = self.h(self.length, 0, s_bb)
        h_s_cc = self.h(self.length, 0, s_cc)
        
        k_val = self.f(1, 0, h_s_bb)
        n_val = self.f(2, 0, h_s_cc)

        oo = raspy_sum(raspy_sum(h_s_aa, k_val), n_val)
        pp = self.i_func(self.length, 1, oo)
        qq = self.j(pp)
        return qq

    def v(self):
        # Addition part
        m_true = self.m("+", True)
        m_false = self.m("+", False)
        
        aa = raspy_sum(self.a(m_true), self.a(m_false))
        
        carry_logic = []
        for x in aa:
            if x > 9: carry_logic.append("1")
            elif x == 9: carry_logic.append("<")
            else: carry_logic.append("0")
            
        bb_pre_n = self.f(-1, "0", carry_logic)
        match_n = [c != "<" for c in bb_pre_n]
        bb_post_n_str = self.n(match_n, bb_pre_n)
        bb_post_n = [int(c) for c in bb_post_n_str]
        
        cc = raspy_mod(raspy_sum(aa, bb_post_n), 10)

        # "pwned" part
        dd = self.u()
        
        ee = [103, 101, 116, 32, 112, 119, 110, 101, 100] + [33] * (self.length - 9)
        
        ff = self.where([d == 1 for d in dd], ee, cc)
        
        cond_aesthetic = [(d == 1 and (i + 1 == self.length) and (i > 10)) for i, d in zip(self.indices, dd)]
        ff = self.where(cond_aesthetic, 49, ff)

        # Convert numbers to ASCII
        ff = self.where([val == 0 for val in ff], 48, ff)
        ff = self.where([val == 1 for val in ff], 49, ff)
        ff = self.where([val == 2 for val in ff], 50, ff)
        ff = self.where([val == 3 for val in ff], 51, ff)
        ff = self.where([val == 4 for val in ff], 52, ff)
        ff = self.where([val == 5 for val in ff], 53, ff)
        ff = self.where([val == 6 for val in ff], 54, ff)
        ff = self.where([val == 7 for val in ff], 55, ff)
        ff = self.where([val == 8 for val in ff], 56, ff)
        ff = self.where([val == 9 for val in ff], 57, ff)

        gg = self.c(ff)
        return gg.lstrip('0') if gg.lstrip('0') else "0"

def solve():
    input1 = "734107+4295754"
    input2 = "5429141+142196"

    runner1 = RaspyRunner(input1)
    output1 = runner1.v()

    runner2 = RaspyRunner(input2)
    output2 = runner2.v()
    
    print(f"{output1};{output2}")

solve()