import textwrap

def explain_solution():
    """
    Analyzes the provided options and explains why one is superior for mitigating
    AitM attacks on DNS in the context of compromised systems.
    """
    best_choice = "G"
    explanation = """
    The challenge describes a scenario where traditional cryptographic defenses against AitM are defeated by adversary implants at the system's core (firmware, kernel). This means the very software running the defenses cannot be trusted.

    Analysis of Solution G: A Secure Blockchain Software Update Mechanism.

    1. Addresses the Root Cause: This solution directly targets how an adversary would install a low-level implant – through a compromised software update or supply chain attack.

    2. Re-establishes Trust: By using a blockchain, every software update can be cryptographically verified against an immutable and tamper-proof ledger. This ensures the integrity of the DNS software itself.

    3. Enables Other Defenses: If the DNS client and server software can be trusted, then the cryptographic protocols they implement (like DNSSEC) can also be trusted to work as intended. This solution fixes the foundation upon which other security measures are built.

    Other options are flawed because they are impractical (B), irrelevant to the protocol (C), protect the wrong component (D, F), are reactive instead of preventative (E), or fail to solve the untrusted system problem (A).
    """

    print(f"The best solution is: {best_choice}")
    print("-" * 20)
    print("\n".join(textwrap.wrap(explanation, width=90)))

explain_solution()