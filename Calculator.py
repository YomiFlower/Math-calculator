import math
import re

try:
    import numpy as np
    has_numpy = True
except ModuleNotFoundError:
    print("Warning: NumPy not found. Using standard math library instead.")
    has_numpy = False

while True:
    # Try to use an interactive backend (TkAgg) if Tkinter is available.
    try:
        import tkinter  # Check if Tkinter is installed.
        import matplotlib
        matplotlib.use('TkAgg')
    except ModuleNotFoundError:
        import matplotlib
        matplotlib.use('Agg')
        print("Warning: Tkinter not found. Using non-interactive backend 'Agg'. Plots will be saved to file instead.")

    import matplotlib.pyplot as plt

    # ======================= HELPER FUNCTIONS =======================
    
    def parse_input(input_str):
        """
        Parse user input that may contain square root expressions.
        Supports formats like:
        - "sqrt(2)" or "√2"
        - "2*sqrt(3)" or "2√3"
        - Regular numbers like "5" or "3.14"
        
        Returns the calculated float value.
        """
        # Remove all spaces
        input_str = input_str.replace(" ", "")
        
        # First, try to parse as a simple float
        try:
            return float(input_str)
        except ValueError:
            pass
        
        # Handle sqrt expressions with Unicode symbol √
        if "√" in input_str:
            # Case 1: Just √n (e.g. "√2")
            if input_str.startswith("√"):
                try:
                    radicand = float(input_str[1:])
                    return math.sqrt(radicand)
                except ValueError:
                    pass
            
            # Case 2: coefficient*√n (e.g. "2√3")
            sqrt_parts = input_str.split("√")
            if len(sqrt_parts) == 2 and sqrt_parts[0] and sqrt_parts[1]:
                try:
                    coefficient = float(sqrt_parts[0])
                    radicand = float(sqrt_parts[1])
                    return coefficient * math.sqrt(radicand)
                except ValueError:
                    pass
        
        # Handle sqrt expressions with sqrt() notation
        sqrt_pattern = r"(\d*)\s*\*?\s*sqrt\(\s*(\d+\.?\d*)\s*\)"
        match = re.match(sqrt_pattern, input_str)
        if match:
            coeff_str, radicand_str = match.groups()
            radicand = float(radicand_str)
            
            # If coefficient is empty, default to 1
            coeff = 1.0
            if coeff_str:
                coeff = float(coeff_str)
                
            return coeff * math.sqrt(radicand)
        
        # If all parsing attempts failed
        raise ValueError(f"Could not parse '{input_str}'. Please use a valid number or square root expression.")

    # ======================= SOLVER FUNCTIONS =======================

    def solve_sss(a, b, c):
        """
        Solve a triangle given all three sides (SSS).
        Sides are: a = BC, b = AC, c = AB.
        Angles: A at vertex A, B at vertex B, C at vertex C.
        """
        if a + b <= c or a + c <= b or b + c <= a:
            raise ValueError("The provided sides do not form a valid triangle.")
        # Calculate angles using the law of cosines
        A = math.acos((b**2 + c**2 - a**2) / (2 * b * c))
        B = math.acos((a**2 + c**2 - b**2) / (2 * a * c))
        C = math.acos((a**2 + b**2 - c**2) / (2 * a * b))
        # Convert angles to degrees
        A_deg = math.degrees(A)
        B_deg = math.degrees(B)
        C_deg = math.degrees(C)
        return a, b, c, A_deg, B_deg, C_deg

    def solve_sas(a, C, b):
        """
        Solve a triangle given two sides and the included angle (SAS).
        a, b are sides and C (in degrees) is the included angle.
        """
        C_rad = math.radians(C)
        c = math.sqrt(a**2 + b**2 - 2 * a * b * math.cos(C_rad))
        A = math.asin(a * math.sin(C_rad) / c)
        B = math.asin(b * math.sin(C_rad) / c)
        return a, b, c, math.degrees(A), math.degrees(B), float(C)

    def solve_asa(A, c, B):
        """
        Solve a triangle given two angles and the included side (ASA).
        c is the included side between angles A and B.
        """
        A_rad = math.radians(A)
        B_rad = math.radians(B)
        C_rad = math.pi - A_rad - B_rad
        a = c * math.sin(A_rad) / math.sin(C_rad)
        b = c * math.sin(B_rad) / math.sin(C_rad)
        return a, b, c, A, B, math.degrees(C_rad)

    def solve_ssa(a, b, A):
        """
        Solve a triangle given two sides and a non-included angle (SSA).
        Returns one solution or two possible solutions.
        """
        A_rad = math.radians(A)
        h = b * math.sin(A_rad)
        if a < h:
            raise ValueError("No valid triangle can be formed with the given sides and angle.")
        elif abs(a - h) < 1e-6:  # right triangle case
            B_rad = math.asin(h / a)
            C_rad = math.pi - A_rad - B_rad
            c = a * math.cos(B_rad)
            return a, b, c, A, math.degrees(B_rad), math.degrees(C_rad)
        else:
            B1_rad = math.asin(b * math.sin(A_rad) / a)
            B2_rad = math.pi - B1_rad
            C1_rad = math.pi - A_rad - B1_rad
            C2_rad = math.pi - A_rad - B2_rad
            c1 = a * math.sin(C1_rad) / math.sin(A_rad)
            c2 = a * math.sin(C2_rad) / math.sin(A_rad)
            A_deg = A
            B1 = math.degrees(B1_rad)
            B2 = math.degrees(B2_rad)
            C1 = math.degrees(C1_rad)
            C2 = math.degrees(C2_rad)
            return [(a, b, c1, A_deg, B1, C1), (a, b, c2, A_deg, B2, C2)]

    def solve_aas(A, B, a):
        """
        Solve a triangle given two angles (A and B) and a non-included side a (opposite angle A).
        Computes the third angle C = 180 - A - B, then uses the law of sines.
        """
        C = 180 - A - B
        if C <= 0:
            raise ValueError("The provided angles do not form a valid triangle.")
        b = a * math.sin(math.radians(B)) / math.sin(math.radians(A))
        c = a * math.sin(math.radians(C)) / math.sin(math.radians(A))
        return a, b, c, A, B, C

    def solve_isosceles_right(side_length, side_type):

        if side_type == 'leg':
            # If we know one of the equal sides (legs)
            a = b = side_length
            c = side_length * math.sqrt(2)  # Hypotenuse = leg * √2
        else:  # side_type == 'hypotenuse'
            # If we know the hypotenuse
            c = side_length
            a = b = side_length / math.sqrt(2)  # Leg = hypotenuse / √2
            
        # The angles in an isosceles right triangle are 45°, 45°, and 90°
        A = B = 45.0
        C = 90.0
        
        return a, b, c, A, B, C

    def solve_angle():
        """
        Find an angle using an inverse trigonometric function (sin, cos, or tan)
        based on a provided ratio.
        """
        print("Select the trigonometric function to compute the angle:")
        print("1. Sine (sin) - angle = arcsin(value)")
        print("2. Cosine (cos) - angle = arccos(value)")
        print("3. Tangent (tan) - angle = arctan(value)")
        func_choice = input("Enter the number corresponding to your choice: ")
        
        if func_choice == '1':
            value = float(input("Enter the ratio value (must be between -1 and 1): "))
            if not (-1 <= value <= 1):
                print("Value out of domain for arcsin.")
                return
            angle = math.degrees(math.asin(value))
            print(f"Angle = {angle:.2f}° (using arcsin)")
        elif func_choice == '2':
            value = float(input("Enter the ratio value (must be between -1 and 1): "))
            if not (-1 <= value <= 1):
                print("Value out of domain for arccos.")
                return
            angle = math.degrees(math.acos(value))
            print(f"Angle = {angle:.2f}° (using arccos)")
        elif func_choice == '3':
            value = float(input("Enter the ratio value: "))
            angle = math.degrees(math.atan(value))
            print(f"Angle = {angle:.2f}° (using arctan)")
        else:
            print("Invalid choice.")

    # ======================= DRAWING FUNCTION =======================

    def draw_triangle_sss(a, b, c, A, B, C):
        """
        Draw the triangle with annotated side lengths and angles.
        Assumes: For SSS (or similar solved triangle), we use:
        - Vertex A at (0,0)
        - Vertex B at (c, 0)  where side AB = c
        - Vertex C computed from side AC = b and angle A at vertex A.
        """
        # Convert angle A (vertex A) to radians
        A_rad = math.radians(A)
        # Coordinates for vertices:
        A_coords = (0, 0)
        B_coords = (c, 0)
        C_coords = (b * math.cos(A_rad), b * math.sin(A_rad))
        
        # Create lists for polygon vertices (closed loop)
        xs = [A_coords[0], B_coords[0], C_coords[0], A_coords[0]]
        ys = [A_coords[1], B_coords[1], C_coords[1], A_coords[1]]
        
        plt.figure()
        plt.plot(xs, ys, 'k-', marker='o')
        plt.title("Triangle")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.axis("equal")
        plt.grid(True)
        
        # Annotate vertices with labels
        plt.text(A_coords[0], A_coords[1], "  A", fontsize=12, color='blue')
        plt.text(B_coords[0], B_coords[1], "  B", fontsize=12, color='blue')
        plt.text(C_coords[0], C_coords[1], "  C", fontsize=12, color='blue')
        
        # Annotate side lengths at midpoints:
        mid_AB = ((A_coords[0] + B_coords[0]) / 2, (A_coords[1] + B_coords[1]) / 2)
        mid_BC = ((B_coords[0] + C_coords[0]) / 2, (B_coords[1] + C_coords[1]) / 2)
        mid_CA = ((C_coords[0] + A_coords[0]) / 2, (C_coords[1] + A_coords[1]) / 2)
        plt.text(*mid_AB, f'{c:.2f}', ha='center', va='bottom', color='red')
        plt.text(*mid_BC, f'{a:.2f}', ha='center', va='bottom', color='red')
        plt.text(*mid_CA, f'{b:.2f}', ha='center', va='top', color='red')
        
        # Annotate angles near each vertex:
        plt.text(A_coords[0], A_coords[1], f'\n{A:.1f}°', ha='left', va='top', color='green')
        plt.text(B_coords[0], B_coords[1], f'\n{B:.1f}°', ha='right', va='top', color='green')
        plt.text(C_coords[0], C_coords[1], f'\n{C:.1f}°', ha='center', va='bottom', color='green')
        
        # Always save to file without asking
        plt.savefig("triangle.png")
        print("Triangle visualization saved as 'triangle.png'.")
        
        # Only show interactive plot if using TkAgg backend
        current_backend = plt.get_backend()
        if current_backend == "TkAgg":
            plt.show()

    # ======================= MAIN FUNCTION =======================

    def display_result(result):
        a, b, c, A, B, C = result
        print(f"\nSides: a = {a:.2f}, b = {b:.2f}, c = {c:.2f}")
        print(f"Angles: A = {A:.2f}°, B = {B:.2f}°, C = {C:.2f}°")

    def main():
        print("Triangle Solver")
        print("Select the known values:")
        print("1. Three sides (SSS)")
        print("2. Two sides and included angle (SAS)")
        print("3. Two angles and included side (ASA)")
        print("4. Two sides and non-included angle (SSA)")
        print("5. Find an angle from a trigonometric ratio")
        print("6. Two angles and a non-included side (AAS)")
        print("7. Isosceles right triangle (one side)")
        
        choice = input("Enter the number corresponding to your choice: ")
        
        try:
            if choice == '1':
                print("Enter side lengths (you can use square root expressions like 'sqrt(2)' or '√2'):")
                a_input = input("Enter side a (BC): ")
                b_input = input("Enter side b (AC): ")
                c_input = input("Enter side c (AB): ")
                
                a = parse_input(a_input)
                b = parse_input(b_input)
                c = parse_input(c_input)
                
                result = solve_sss(a, b, c)
                display_result(result)
                draw_triangle_sss(a, b, c, result[3], result[4], result[5])
            elif choice == '2':
                print("Enter values (you can use square root expressions for sides like 'sqrt(2)' or '√2'):")
                a_input = input("Enter side a: ")
                C = float(input("Enter included angle C (in degrees): "))
                b_input = input("Enter side b: ")
                
                a = parse_input(a_input)
                b = parse_input(b_input)
                
                result = solve_sas(a, C, b)
                display_result(result)
                draw_triangle_sss(result[0], result[1], result[2], result[3], result[4], result[5])
            elif choice == '3':
                print("Enter values (you can use square root expressions for sides like 'sqrt(2)' or '√2'):")
                A = float(input("Enter angle A (in degrees): "))
                c_input = input("Enter included side c: ")
                B = float(input("Enter angle B (in degrees): "))
                
                c = parse_input(c_input)
                
                result = solve_asa(A, c, B)
                display_result(result)
                draw_triangle_sss(result[0], result[1], result[2], result[3], result[4], result[5])
            elif choice == '4':
                print("Enter values (you can use square root expressions for sides like 'sqrt(2)' or '√2'):")
                a_input = input("Enter side a: ")
                b_input = input("Enter side b: ")
                A = float(input("Enter angle A (in degrees): "))
                
                a = parse_input(a_input)
                b = parse_input(b_input)
                
                result = solve_ssa(a, b, A)
                if isinstance(result, list):
                    print("\nTwo possible solutions exist:")
                    for i, res in enumerate(result, 1):
                        print(f"\nSolution {i}:")
                        display_result(res)
                        draw_triangle_sss(res[0], res[1], res[2], res[3], res[4], res[5])
                        plt.savefig(f"triangle_solution_{i}.png")
                        print(f"Triangle visualization for solution {i} saved as 'triangle_solution_{i}.png'.")
                else:
                    display_result(result)
                    draw_triangle_sss(result[0], result[1], result[2], result[3], result[4], result[5])
            elif choice == '5':
                solve_angle()
            elif choice == '6':
                print("Solving triangle using AAS (Two angles and a non-included side).")
                print("Enter values (you can use square root expressions for sides like 'sqrt(2)' or '√2'):")
                A = float(input("Enter angle A (in degrees, where side a is opposite A): "))
                B = float(input("Enter angle B (in degrees): "))
                a_input = input("Enter side a (opposite angle A): ")
                
                a = parse_input(a_input)
                
                result = solve_aas(A, B, a)
                display_result(result)
                draw_triangle_sss(result[0], result[1], result[2], result[3], result[4], result[5])
            elif choice == '7':
                print("Solving isosceles right triangle.")
                print("In an isosceles right triangle, two sides are equal and one angle is 90°.")
                print("You can use square root expressions like 'sqrt(2)' or '√2'.")
                side_type = input("Do you know the length of a leg or the hypotenuse? (leg/hypotenuse): ").strip().lower()
                if side_type not in ['leg', 'hypotenuse']:
                    print("Invalid input. Please enter 'leg' or 'hypotenuse'.")
                    return
                side_input = input(f"Enter the length of the {side_type}: ")
                
                side_length = parse_input(side_input)
                
                result = solve_isosceles_right(side_length, side_type)
                display_result(result)
                draw_triangle_sss(result[0], result[1], result[2], result[3], result[4], result[5])
            else:
                print("Invalid choice.")
        except ValueError as e:
            print(f"Error: {e}")

    if __name__ == "__main__":
        main()
