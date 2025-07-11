import socket

def simulate_multi_modal_dns_verification():
    """
    This script simulates the Multi-Modal Verification Process to detect a
    potential Adversary-in-the-Middle (AitM) attack on DNS.
    """
    domain = "example.com"
    # We will use the real IP for the 'trusted' paths
    try:
        # Get a real IP address for our simulation of trusted paths
        real_ip = socket.gethostbyname(domain)
    except socket.gaierror:
        # Fallback if there is no internet connection
        real_ip = "93.184.216.34"

    # Attacker's IP for the 'compromised' path
    fake_ip = "10.13.37.10"

    # Simulate making DNS queries over three independent communication paths.
    # The primary path is compromised in this scenario.
    paths = {
        "Path 1 (Primary Wired Network)": fake_ip,
        "Path 2 (Guest Wi-Fi Network)": real_ip,
        "Path 3 (Cellular 5G Network)": real_ip
    }

    print(f"--- Starting Multi-Modal DNS Verification for '{domain}' ---")
    
    resolved_ips = {}
    for path_name, ip in paths.items():
        print(f"Querying via {path_name}... resolved to IP: {ip}")
        resolved_ips[path_name] = ip
        
    print("\n--- Cross-Validating Results ---")

    ip_values = list(resolved_ips.values())
    
    # This check represents our "final equation" to verify integrity.
    # The equation checks if all resolved IP addresses are equal.
    print("Verification Equation:")
    # We build and print the equation string with the actual numbers (IPs).
    equation_string = " == ".join(ip_values)
    print(f"  {equation_string}")

    # Evaluate the condition
    is_consistent = len(set(ip_values)) == 1

    print(f"Result of Equation: {is_consistent}\n")

    if is_consistent:
        print("✅ SUCCESS: DNS resolution is consistent across all paths.")
        print(f"Verified IP for '{domain}' is {ip_values[0]}")
    else:
        print("🚨 WARNING: Inconsistent DNS results detected!")
        print("This indicates a high probability of an AitM attack.")
        for path_name, ip in resolved_ips.items():
             print(f"  - {path_name}: {ip}")

if __name__ == "__main__":
    simulate_multi_modal_dns_verification()
