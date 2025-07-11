def get_kappa_definition():
    """
    This function provides the definition of the constant kappa based on the
    analysis of the provided synaptic plasticity model.
    """
    # The definition of kappa is derived from the stability analysis of the
    # synaptic weight equation. It is a dimensionless ratio of terms related
    # to the presynaptic accumulator v_k.

    # Symbolic representation of the statistical moments of v_k and the offset rho
    mean_v_k = "<v_k>"
    variance_v_k = "Var(v_k)"
    rho = "rho"

    # The numerator represents the mean-driven component of plasticity for a single synapse
    numerator = f"-{mean_v_k}*({mean_v_k} + {rho})"

    # The denominator represents the fluctuation-driven component
    denominator = variance_v_k

    # Print the full definition
    print("The definition of kappa is given by the following formula:")
    print(f"kappa = ({numerator}) / ({denominator})")
    print("\nWhere:")
    print(f"  {mean_v_k}: The mean of the presynaptic accumulator v_k.")
    print(f"  {variance_v_k}: The variance of the presynaptic accumulator v_k.")
    print(f"  {rho}: The offset constant in the weight update rule.")

if __name__ == "__main__":
    get_kappa_definition()