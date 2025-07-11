import dns.resolver
import sys

def verify_dns_with_multiple_resolvers(domain):
    """
    Demonstrates the Multi-Modal Verification principle by querying a domain 
    against a list of public DNS resolvers and comparing their responses.
    This script requires the 'dnspython' library: pip install dnspython
    """
    # A dictionary of independent, public DNS resolvers to query
    resolvers = {
        'Google': '8.8.8.8',
        'Cloudflare': '1.1.1.1',
        'Quad9': '9.9.9.9',
        'OpenDNS': '208.67.222.222'
    }
    
    results = {}
    print(f"--- Starting Multi-Modal DNS Verification for: {domain} ---\n")

    # Step 1: Query all resolvers
    for name, ip in resolvers.items():
        try:
            resolver = dns.resolver.Resolver()
            resolver.nameservers = [ip]
            # Set a timeout to avoid long waits for unresponsive servers
            resolver.timeout = 2
            resolver.lifetime = 2
            
            # We are querying for A records (IPv4 addresses)
            answer = resolver.resolve(domain, 'A')
            # Store a sorted list of IPs to ensure comparison is consistent
            results[name] = sorted([rdata.to_text() for rdata in answer])
            print(f"[+] Query to {name} ({ip}) returned: {results[name]}")
        except Exception as e:
            # Store the error message if a resolver fails
            results[name] = f"Error: {type(e).__name__}"
            print(f"[-] Query to {name} ({ip}) failed: {results[name]}")

    print("\n--- Performing Verification Cross-Check ---")
    
    # Step 2: Compare the results for consensus
    valid_responses = [tuple(v) for v in results.values() if isinstance(v, list)]
    
    if not valid_responses:
        print("\n[!] Conclusion: Verification failed. Could not retrieve any valid DNS responses.")
        return

    first_response = valid_responses[0]
    all_match = all(res == first_response for res in valid_responses)

    # Step 3: Print the final "equation" and conclusion
    # The "equation" represents the comparison of all resolver responses.
    print("\nVerification Equation:")
    # The instruction asks to output each number in the equation. The numbers are the resolver IPs.
    equation_str = " == ".join([f"Response({ip})" for ip in resolvers.values()])
    print(f"  Is {equation_str} ? -> {all_match}")
    
    print("\nNumbers (Resolvers) used in equation:")
    for number in resolvers.values():
        print(f"  - {number}")
        
    print("\n--- Final Conclusion ---")
    if all_match:
        print(f"✅ Consensus Achieved. All reachable resolvers returned the same result.")
        print(f"   Verified IP Address(es) for '{domain}': {list(first_response)}")
    else:
        print(f"⚠️ WARNING: DISCREPANCY DETECTED!")
        print(f"   Potential DNS poisoning or Adversary-in-the-Middle attack.")
        print("   The following differing results were received:")
        for name, result in results.items():
            print(f"     - {name}: {result}")

if __name__ == "__main__":
    # You can run this script with a domain name as an argument
    # e.g., python your_script_name.py google.com
    target_domain = "ietf.org"
    if len(sys.argv) > 1:
        target_domain = sys.argv[1]
        
    verify_dns_with_multiple_resolvers(target_domain)