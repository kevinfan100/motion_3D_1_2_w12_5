---
description: Run lambda_c sweep on a verification script
argument: script name (e.g., verify_eq13)
---

Run the specified verification script using MATLAB MCP and format the results as a table showing variance and recovered a_x for each lambda_c value.

If the script name is not provided, ask the user which script to run.

Present results in a formatted table:

| lambda_c | C(lc) | sigma^2 (um^2) | a_xm (um/pN) | Error (%) |
