import shutil
import os
from pathlib import Path
 
def migrate_database():
    old_db_path = "gap_fade_trades.db"
    new_db_dir = Path("/data")
    new_db_path = new_db_dir / "gap_fade_trades.db"
    
    new_db_dir.mkdir(parents=True, exist_ok=True)
    
    if os.path.exists(old_db_path):
        if not os.path.exists(str(new_db_path)):
            try:
                print(f"[MIGRATION] Found {old_db_path}")
                print(f"[MIGRATION] Copying to {new_db_path}...")
                shutil.copy2(old_db_path, str(new_db_path))
                print(f"[MIGRATION] ✓ Success: {old_db_path} → {new_db_path}")
            except Exception as e:
                print(f"[MIGRATION] ✗ Error: {e}")
                raise
        else:
            print(f"[MIGRATION] {new_db_path} already exists, skipping")
    else:
        print(f"[MIGRATION] No {old_db_path} found (first run?), continuing...")
 
if __name__ == "__main__":
    migrate_database()
