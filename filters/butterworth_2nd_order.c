#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846

double butterworth_2nd_order(double x, double y[], double dt, double cutoff) {
    static double a[3];
    static double b[3];
    static double x_prev[3] = {0};
    static double y_prev[3] = {0};
    
    double wc = 2.0 * PI * cutoff;
    double n = 2.0 * dt * wc;
    double d = n * n + sqrt(2.0) * n + 1;
    
    a[0] = n * n / d;
    a[1] = 2.0 * a[0];
    a[2] = a[0];
    b[1] = 2.0 * (n * n - 1) / d;
    b[2] = (n * n - sqrt(2.0) * n + 1) / d;
    
    y[0] = a[0] * x + a[1] * x_prev[1] + a[2] * x_prev[0] - b[1] * y_prev[1] - b[2] * y_prev[0];
    
    x_prev[0] = x_prev[1];
    x_prev[1] = x;
    y_prev[0] = y_prev[1];
    y_prev[1] = y[0];
    
    return y[0];
}

int main() {
    double data[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double filtered_data[10];
    double dt = 0.1;
    double cutoff = 1.0;
    
    for (int i = 0; i < 10; i++) {
        filtered_data[i] = butterworth_2nd_order(data[i], &filtered_data[i], dt, cutoff);
        printf("%.6f\n", filtered_data[i]);
    }
    
    return 0;
}
