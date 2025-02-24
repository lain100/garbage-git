public class foo {
  public static void main(String[] args) {
    double x0 = 0;
    double[] y0 = {0, 1};
    showAnswer( rungeKutta(x0, y0, 0.5, 200) );
  }
  
  public static double[][] rungeKutta(double x0, double y0, double h, int s) {
    double[] x = new double[s + 1];
    double[] y = new double[s + 1];
    x[0] = x0;
    y[0] = y0;
    
    for (int i = 0; i < s; i++) {
      double k1 = h * f(x[i], y[i]);
      double k2 = h * f(x[i] + h/2.0, y[i] + k1/2.0);
      double k3 = h * f(x[i] + h/2.0, y[i] + k2/2.0);
      double k4 = h * f(x[i] + h, y[i] + k3);
      x[i + 1] = x[i] + h;
      y[i + 1] = y[i] + (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;
    }
    
    double[][] ans = new double[s + 1][2];
    for (int i = 0; i < s + 1; i++) {
      ans[i][0] = x[i];
      ans[i][1] = y[i];
    }
    
    return ans;
  }
  
  public static double[][] rungeKutta(double x0, double[] y0, double h, int s) {
    int n = y0.length;
    double[] x = new double[s + 1];
    double[][] y = new double[s + 1][n];
    double[][] k = new double[n][4];
    double[] tmp = new double[n];
    double[] ratio = {0.5, 0.5, 1};
    
    x[0] = x0;
    for (int j = 0; j < n; j++) y[0][j] = y0[j];
    
    for (int i = 0; i < s; i++) {
      x[i + 1] = x[i] + h;
      
      double[] f = f(x[i], y[i]);
      for (int j = 0; j < n; j++) k[j][0] = h * f[j];
      
      for (int m = 0; m < 3; m++) {
        for (int j = 0; j < n; j++) tmp[j] = y[i][j] + ratio[m] * k[j][m];
        f = f(x[i] + ratio[m] * h, tmp);
        for (int j = 0; j < n; j++) k[j][m + 1] = h * f[j];
      }
      for (int j = 0; j < n; j++) {
        y[i + 1][j] = y[i][j] + (k[j][0] + 2*k[j][1] + 2*k[j][2] + k[j][3])/6.0;
      }
    }
    
    double[][] ans = new double[s + 1][2];
    for (int i = 0; i < s + 1; i++) {
      ans[i][0] = x[i];
      ans[i][1] = y[i][0];
    }
    
    return ans;
  }
  
  public static double[][] heun(double x0, double y0, double h, int s) {
    double[] x = new double[s + 1]; 
    double[] y = new double[s + 1]; 
    x[0] = x0;
    y[0] = y0;
    
    for (int i = 0; i < s; i++) {
      double k1 = h * f(x[i], y[i]);
      double k2 = h * f(x[i] + h, y[i] + k1);
      x[i + 1] = x[i] + h;
      y[i + 1] = y[i] + (k1 + k2) / 2.0;
    }
    
    double[][] ans = new double[s + 1][2];
    for (int i = 0; i < s + 1; i++) {
      ans[i][0] = x[i];
      ans[i][1] = y[i];
    }
    
    return ans;
  }
  
  public static double[][] heun(double x0, double[] y0 ,double h, int s) {
    int n = y0.length;
    double[] x = new double[s + 1];
    double[][] y = new double[s + 1][n];
    double[][] k = new double[n][2];
    double[] tmp = new double[n];
    
    x[0] = x0;
    for (int j = 0; j < n; j++) y[0][j] = y0[j];
    
    for (int i = 0; i < s; i++) {
      x[i + 1] = x[i] + h;
      
      double[] f = f(x[i], y[i]);
      for (int j = 0; j < n; j++) k[j][0] = h * f[j];
      
      for (int j = 0; j < n; j++) tmp[j] = y[i][j] + k[j][0];
      f = f(x[i] + h, tmp);
      for (int j = 0; j < n; j++) k[j][1] = h * f[j];
      
      for (int j = 0; j < n; j++) {
        y[i + 1][j] = y[i][j] + (k[j][0] + k[j][1]) / 2.0;
      }
    }
    
    double[][] ans = new double[s + 1][2];
    for (int i = 0; i < s + 1; i++) {
      ans[i][0] = x[i];
      ans[i][1] = y[i][0];
    }
    
    return ans;
  }
  
  public static double[][] euler(double x0, double y0, double h, int s) {
    double[][] ans = new double[s + 1][2];
    ans[0][0] = x0;
    ans[0][1] = y0;
    
    for (int i = 0; i < s; i++) {
      ans[i + 1][0] = ans[i][0] + h;
      ans[i + 1][1] = ans[i][1] + h*f(ans[i][0], ans[i][1]);
    }
    
    return ans;
  }
  
  public static double[][] euler(double x0, double[] y0, double h, int s) {
    int n = y0.length;
    double[] x = new double[s + 1];
    double[][] y = new double[s + 1][n];
    
    x[0] = x0;
    for (int j = 0; j < n; j++) y[0][j] = y0[j];
    
    for (int i = 0; i < s; i++) {
      double[] f = f(x[i], y[i]);

      x[i + 1] = x[i] + h;
      for (int j = 0; j < n; j++) {
        y[i + 1][j] = y[i][j] + h*f[j];
      }
    }
    
    double [][] ans = new double[s + 1][2];
    for (int i = 0; i < s + 1; i++) {
      ans[i][0] = x[i];
      ans[i][1] = y[i][0];
    }
    
    return ans;
  }
  
  public static void showAnswer(double[][] ans) {
    int m = ans.length;
    int n = ans[0].length;
    
    for (int i = 0; i < m; i++) {      
      for (int j = 0; j < n; j++) {
        System.out.printf("%.8f ", ans[i][j]);
      }
      System.out.println();
    }
  }
  
  public static double f(double x, double y) {
    return 3*x*x*y;
  }
  
  public static double[] f(double x, double[] y) {
    double[] ans = new double[3];
    
    ans[0] = y[1];                    // y0'(x) = y1(x)
    ans[1] = -1.01*y[0] -0.2*y[1];    // y1'(x) = f(x, y)
    
    return ans;
  }
  
}
