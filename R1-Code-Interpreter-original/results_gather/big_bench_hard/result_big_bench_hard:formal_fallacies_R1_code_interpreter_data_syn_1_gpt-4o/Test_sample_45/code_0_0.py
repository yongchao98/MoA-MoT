def is_valid():
    # Define the premise and conclusion
    # Premise: If not F, then not N
    # Conclusion: If N, then F
    # We need to check if the contrapositive of the premise is equivalent to the conclusion

    # Logical equivalence check
    # If not F then not N is equivalent to If N then F
    # This is a tautology check: (not F -> not N) is equivalent to (N -> F)
    
    # Let's assume F and N are boolean variables
    # We will check all possible truth values for F and N
    for F in [True, False]:
        for N in [True, False]:
            # Premise: If not F, then not N
            premise = (not F) <= (not N)
            # Conclusion: If N, then F
            conclusion = N <= F
            # Check if they are equivalent
            if premise != conclusion:
                return False
    return True

# Check if the argument is valid
print("valid" if is_valid() else "invalid")