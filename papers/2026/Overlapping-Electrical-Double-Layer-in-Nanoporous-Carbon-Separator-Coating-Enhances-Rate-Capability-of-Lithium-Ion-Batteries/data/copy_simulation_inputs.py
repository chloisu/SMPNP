import argparse
import pathlib
import shutil
import sys

def copy_simulation_inputs(source_dir_str, target_dir_str):
    """
    Copies 'input.yml' and 'parameters.yml' from source to a target directory.

    SAFETY: Aborts if the target directory already exists and is not empty
    to prevent accidental data loss.
    """
    
    source_dir = pathlib.Path(source_dir_str)
    target_dir = pathlib.Path(target_dir_str)
    
    files_to_copy = ["input.yml", "parameters.yml"]

    # --- CRITICAL SAFETY CHECK ---
    if target_dir.exists() and any(target_dir.iterdir()):
        print(f"⛔️ ERROR: Destination directory '{target_dir}' already exists and is not empty.")
        print("To prevent accidental data loss, this script will not modify an existing directory.")
        print("\nACTION REQUIRED:")
        print(f"  Please manually delete or move the '{target_dir}' folder and run the script again.")
        sys.exit(1)

    # --- Input Validation and Setup ---
    if not source_dir.is_dir():
        print(f"Error: Source directory '{source_dir}' not found or is not a directory.")
        sys.exit(1)

    target_dir.mkdir(parents=True, exist_ok=True)
    print(f"Source:      '{source_dir}'")
    print(f"Destination: '{target_dir}' (Verified as safe to write to)")
    print("-" * 50)

    # --- Main Logic: Iterate and Copy ---
    case_directories = [d for d in source_dir.iterdir() if d.is_dir()]

    if not case_directories:
        print("Warning: No case subdirectories found in the source directory.")
        return

    incomplete_cases = {}
    copied_file_count = 0
    
    for case_dir in sorted(case_directories):
        case_name = case_dir.name
        print(f"Processing case: {case_name}")

        found_files_for_this_case = []
        missing_files_for_this_case = []

        for filename in files_to_copy:
            if (case_dir / filename).exists():
                found_files_for_this_case.append(filename)
            else:
                missing_files_for_this_case.append(filename)
        
        if found_files_for_this_case:
            dest_case_dir = target_dir / case_name
            dest_case_dir.mkdir(exist_ok=True)

            for filename in found_files_for_this_case:
                shutil.copy2(case_dir / filename, dest_case_dir / filename)
                print(f"  -> Copied '{filename}'")
                copied_file_count += 1
            
            if missing_files_for_this_case:
                incomplete_cases[case_name] = missing_files_for_this_case
                for filename in missing_files_for_this_case:
                    print(f"  -> Warning: '{filename}' not found.")
        else:
            print("  -> Skipped: No target files found.")

    # --- Final Summary ---
    print("-" * 50)
    print(f"Script finished. Copied a total of {copied_file_count} file(s).")

    if incomplete_cases:
        print("\n--- WARNING SUMMARY ---")
        print("The following cases were incomplete (but files found were still copied):")
        for case, missing_list in incomplete_cases.items():
            # --- THIS IS THE CORRECTED PART ---
            # Create the formatted string of missing files separately.
            missing_files_str = ", ".join(f"'{f}'" for f in missing_list)
            # Now use the simple string in the final print statement.
            print(f"  - {case}: Was missing {missing_files_str}")
    else:
        print("\nAll cases contained all required input files or had no files to copy.")
    
    sys.exit(0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Safely copy simulation input files to a new destination directory.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument("source", help="Path to the source directory containing case subfolders.")
    parser.add_argument("destination", help="Path to the new destination directory to be created.")

    args = parser.parse_args()
    copy_simulation_inputs(args.source, args.destination)