import math

# Code below written and checked with help of ChatGPT
def t_distribution_probability(degrees_freedom, z_value):
    """
    Calculate the probability of the t-distribution given degrees of freedom and z value.

    Args:
        degrees_freedom (int): Degrees of freedom for the t-distribution.
        z_value (float): Value of z for which the probability is calculated.

    Returns:
        float: Probability of the t-distribution.
    """
    # Calculate the t-distribution probability using the formula
    numerator = math.gamma((degrees_freedom + 1) / 2)
    denominator = math.sqrt(degrees_freedom * math.pi) * math.gamma(degrees_freedom / 2)
    coefficient = math.pow(1 + (z_value ** 2 / degrees_freedom), -(degrees_freedom + 1) / 2)
    return (numerator / denominator) * coefficient


def main():
    """Main function to get user input and calculate t-distribution probability."""
    print("Enter degrees of freedom (integer):")
    degrees_freedom = int(input())

    print("Enter z value (float):")
    z_value = float(input())

    probability = t_distribution_probability(degrees_freedom, z_value)

    print("Probability:", probability)


if __name__ == "__main__":
    main()
