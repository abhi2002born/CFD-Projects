#include <bits/stdc++.h>
#include <chrono>

using namespace std;
using namespace std::chrono;

#define E 1e-5 // Convergence tolerance
#define m 128  // Grid size (128x128)



void Matrix_A_Vect_Multi(const vector<double> &x, vector<double> &result)
{
    int nx = m;
    result.assign(nx * nx, 0.0);
    
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            int idx = i * nx + j;
            double val = 4.0 * x[idx];
            
            if (j > 0)    val -= x[idx - 1]; // West neighbor
            if (j < nx-1) val -= x[idx + 1]; // East neighbor
            if (i > 0)    val -= x[idx - nx]; // South neighbor
            if (i < nx-1) val -= x[idx + nx]; // North neighbor
            
            // Adjust boundary conditions
            if (i == 0 || i == nx - 1) val += x[idx];
            if (j == 0 || j == nx - 1) val += x[idx];
            
            result[idx] = val;
        }
    }
}

// Function to compute dot product of two vectors
double Dot_Product(const vector<double> &a, const vector<double> &b)
{
    double result = 0.0;
    for (size_t i = 0; i < a.size(); i++)
    {
        result += a[i] * b[i];
    }
    return result;
}

// Function to compute 2-norm of a vector
double norm_2(const vector<double> &a)
{
    double result = 0.0;
    for (size_t i = 0; i < a.size(); i++)
    {
        result += a[i] * a[i];
    }
    return sqrt(result);
}

// Function to compute analytical solution
double analytical_solution(double x, double y, int terms)
{
    double T = 0.0;
    for (int n = 1; n <= terms; n++)
    {
        double coef = (2.0 / M_PI) * ((pow(-1.0, n + 1) + 1.0) / n);
        T += coef * sin(n * M_PI * x) * sinh(n * M_PI * y) / sinh(n * M_PI);
    }
    return T;
}

// Function to save results to file
void save_results(const vector<double> &T, const string &filename)
{
    ofstream file(filename);
    file << "VARIABLES = X, Y, T\n";
    file << "ZONE I=" << m << ", J=" << m << ", F=POINT\n";

    double dx = 1.0 / m;
    double dy = 1.0 / m;

    for (int j = 0; j < m; j++)
    {
        double y = (j + 0.5) * dy;
        for (int i = 0; i < m; i++)
        {
            double x = (i + 0.5) * dx;
            int idx = i + j * m;
            file << fixed << setprecision(9) << x << "\t" << y << "\t" << T[idx] << "\n";
        }
    }
    file.close();
}

// Function to save temperature along centerlines
void save_centerlines(const vector<double> &T, const string &h_filename, const string &v_filename)
{
    ofstream h_file(h_filename);
    ofstream v_file(v_filename);

    double dx = 1.0 / m;
    double dy = 1.0 / m;

    h_file << "VARIABLES = Y, T\n";
    h_file << "ZONE I=" << m << ", F=POINT\n";

    v_file << "VARIABLES = X, T\n";
    v_file << "ZONE I=" << m << ", F=POINT\n";

    // Horizontal centerline (x = 0.5)
    for (int j = 0; j < m; j++)
    {
        double y = (j + 0.5) * dy;
        int idx = m / 2 + j * m;
        h_file << fixed << setprecision(9) << y << "\t" << T[idx] << "\n";
    }

    // Vertical centerline (y = 0.5)
    for (int i = 0; i < m; i++)
    {
        double x = (i + 0.5) * dx;
        int idx = i + (m / 2) * m;
        v_file << fixed << setprecision(9) << x << "\t" << T[idx] << "\n";
    }

    h_file.close();
    v_file.close();
}

int main()
{   
    auto start = high_resolution_clock::now();  // Start timer
    // Problem size
    int n = m * m;

    // Vectors for CG method
    vector<double> r(n, 0.0), d(n, 0.0), b(n, 0.0), temp(n, 0.0), x(n, 0.0);
    vector<double> residual_history;

    // Set up right-hand side vector b with boundary conditions
    // T=0 at left, right, and bottom boundaries
    // T=1 at top boundary
    double dx = 1.0 / m;
    double dy = 1.0 / m;

    // Initialize b with boundary conditions
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            int idx = i + j * m;

            // Interior points
            b[idx] = 0.0;

            // Apply boundary conditions to b
            // Top boundary (T=1)
            if (j == m - 1)
            {
                b[idx] += 2.0; // Contribution from top boundary
            }
        }
    }

    // Initial guess for x is all zeros
    for (int i = 0; i < n; i++)
    {
        x[i] = 0.0;
    }

    // Initialize residual r = b - Ax
    Matrix_A_Vect_Multi(x, temp);
    for (int i = 0; i < n; i++)
    {
        r[i] = b[i] - temp[i];
        d[i] = r[i]; // Initial search direction
    }

    double rs_old = Dot_Product(r, r);
    double rn_old = norm_2(r);
    residual_history.push_back(rn_old);

    // CG iteration
    int iter = 0;
    int max_iter = 10000;
    double alpha = 0.0, beta = 0.0, rs_new = 0.0, rn_new = 0.0;

    cout << "Starting Conjugate Gradient iterations..." << endl;
    cout << "Iteration\tResidual" << endl;
    cout << iter << "\t\t" << rn_old << endl;

    while (rn_old > E && iter < max_iter)
    {
        iter++;

        // Compute alpha
        Matrix_A_Vect_Multi(d, temp);
        alpha = rs_old / Dot_Product(d, temp);

        // Update solution
        for (int i = 0; i < n; i++)
        {
            x[i] += alpha * d[i];
        }

        // Update residual
        if (iter % 50 == 0)
        {
            // Recompute residual periodically to avoid round-off error accumulation
            Matrix_A_Vect_Multi(x, temp);
            for (int i = 0; i < n; i++)
            {
                r[i] = b[i] - temp[i];
            }
        }
        else
        {
            // Standard residual update
            for (int i = 0; i < n; i++)
            {
                r[i] -= alpha * temp[i];
            }
        }

        // Compute new residual norm
        rs_new = Dot_Product(r, r);
        rn_new = norm_2(r);
        residual_history.push_back(rn_new);

        // Compute beta
        beta = rs_new / rs_old;
        rs_old = rs_new;

        // Update search direction
        for (int i = 0; i < n; i++)
        {
            d[i] = r[i] + beta * d[i];
        }

        rn_old = rn_new;

        if (iter % 100 == 0 || rn_new <= E)
        {
            cout << iter << "\t\t" << rn_new << endl;
        }
    }

    cout << "CG converged in " << iter << " iterations." << endl;
    cout << "Final residual: " << rn_new << endl;

    // Save results
    save_results(x, "temperature_cg.dat");
    save_centerlines(x, "T_centerline_H.dat", "T_centerline_V.dat");

    // Compute analytical solution for comparison
    vector<double> T_analytical(n, 0.0);
    for (int j = 0; j < m; j++)
    {
        double y = (j + 0.5) * dy;
        for (int i = 0; i < m; i++)
        {
            double x_pos = (i + 0.5) * dx;
            int idx = i + j * m;
            T_analytical[idx] = analytical_solution(x_pos, y, 100); // Use 100 terms for summation
        }
    }

    // Save analytical solution and print centreline
    save_results(T_analytical, "temperature_analytical.dat");

    save_centerlines(T_analytical, "TA_centerline_H.dat", "TA_centerline_V.dat");

    // Compute error between numerical and analytical solutions
    double max_error = 0.0;
    int x_index, y_index;
    for (int i = 0; i < n; i++)
    {
        double error = fabs(x[i] - T_analytical[i]);
        if (error > max_error)
        {
            max_error = error;
            x_index = i / m;
            y_index = i % m;
        }
    }

    cout << "Maximum error between numerical and analytical solutions: " << max_error << " at i =" << x_index << " and j =" << y_index << endl;

    // Save residual history
    ofstream res_file("residual_history.dat");
    res_file << "VARIABLES = Iteration, Residual\n";
    res_file << "ZONE I=" << residual_history.size() << ", F=POINT\n";
    for (size_t i = 0; i < residual_history.size(); i++)
    {
        res_file << i << "\t" << residual_history[i] << "\n";
    }
    res_file.close();

    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "Execution time: " << duration.count() << " ms" << endl;

    return 0;
}
