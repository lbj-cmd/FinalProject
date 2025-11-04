import sqlite3

# Connect to the database
conn = sqlite3.connect('kirkwood_gap_simulations.db')
cursor = conn.cursor()

# Get list of tables
cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
tables = cursor.fetchall()

print("数据库表:")
for table in tables:
    print(f"  - {table[0]}")

# Get table schemas
for table in tables:
    table_name = table[0]
    print(f"\n{table_name}表结构:")
    cursor.execute(f"PRAGMA table_info({table_name})")
    columns = cursor.fetchall()
    for column in columns:
        print(f"  - {column[1]} ({column[2]})")

# Close connection
conn.close()