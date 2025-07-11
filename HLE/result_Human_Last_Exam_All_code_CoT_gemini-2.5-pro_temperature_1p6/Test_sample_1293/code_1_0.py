import dns.resolver
from collections import Counter

def multi_modal_dns_check(domain_to_check):
    """
    Demonstrates a multi-modal DNS verification process to detect discrepancies
    that could indicate an Adversary-in-the-Middle (AitM) attack.
    """
    print(f"--- Starting Multi-Modal DNS Check for '{domain_to_check}' ---")
    
    # Simulate making requests through different devices/networks/resolvers
    # We include a simulated "compromised" path that returns a malicious IP.
    resolvers = {
        "Google_DNS_on_Primary_Network": "8.8.8.8",
        "Cloudflare_DNS_on_Secondary_Network": "1.1.1.1",
        "Quad9_DNS_on_VPN": "9.9.9.9",
        "Local_ISP_DNS_or_Compromised_Path": "SIMULATED_COMPROMISE" # This will be our fake response
    }
    
    malicious_ip = "10.66.6.1" # A typical non-routable IP used for malicious redirection
    results = {}
    ip_list = []

    print("\nStep 1: Querying multiple DNS resolvers simultaneously...")
    for name, ip in resolvers.items():
        try:
            current_ips = set()
            print(f"  Querying via '{name}'...")
            
            if ip == "SIMULATED_COMPROMISE":
                # This simulates a compromised system/router feeding a malicious IP
                current_ips.add(malicious_ip)
                print(f"    -> Received SUSPICIOUS response: {malicious_ip}")
            else:
                resolver = dns.resolver.Resolver()
                resolver.nameservers = [ip]
                answer = resolver.resolve(domain_to_check, 'A')
                for record in answer:
                    current_ips.add(record.to_text())
                print(f"    -> Received response: {', '.join(current_ips)}")

            results[name] = sorted(list(current_ips))
            ip_list.extend(results[name])
        except Exception as e:
            print(f"    -> Query failed for {name}: {e}")
            results[name] = ["Error"]

    print("\nStep 2: Cross-validating all responses...")
    if not ip_list:
        print("Could not retrieve any valid IP addresses.")
        return

    # Use a counter to find the most common response
    ip_counts = Counter(tuple(ips) for ips in results.values() if ips != ["Error"])
    
    if not ip_counts:
        print("Cross-validation failed. No successful responses to compare.")
        return

    consensus_ips, consensus_count = ip_counts.most_common(1)[0]
    
    print(f"\n--- Verification Analysis ---")
    print(f"Consensus (most common) IP address(es) are: {', '.join(consensus_ips)} ({consensus_count} times)")

    print("\nIndividual Path Results:")
    all_clear = True
    for name, ips in results.items():
        if tuple(ips) == consensus_ips:
            print(f"  [OK] '{name}' returned the consensus IP(s).")
        elif ips == ["Error"]:
            print(f"  [FAIL] '{name}' failed to respond.")
            all_clear = False
        else:
            print(f"  [ALERT!] '{name}' returned a different IP: {', '.join(ips)}. This path may be compromised!")
            all_clear = False
            
    if all_clear:
        print("\nConclusion: All paths returned the same address. The connection appears secure.")
    else:
        print(f"\nConclusion: Discrepancy detected! The IP '{malicious_ip}' is highly suspicious. An AitM attack is possible.")


if __name__ == "__main__":
    # The domain we want to verify.
    target_domain = "www.google.com"
    multi_modal_dns_check(target_domain)