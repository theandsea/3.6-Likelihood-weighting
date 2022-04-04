import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.UUID;

public class likelihood {
    public static void main(String[] args) throws IOException {
        int[][][] graph = graph_new(1000, 700, new int[]{255, 255, 255});
        int times = 8;
        int l = 20;//26;
        int[] x = new int[l];
        double[][] y = new double[4][l];
        double[] thisestimate = null;
        for (int i = 0; i < l; i++, times *= 2) {
            x[i] = times;
            thisestimate = likelihood_s(times);
            y[0][i] = thisestimate[1];
            y[1][i] = thisestimate[4];
            y[2][i] = thisestimate[7];
            y[3][i] = thisestimate[9];
        }
        // x-axis
        //graph_line_relative(graph,0,0,x[l-1],0,3,new int[]{0,0,0},-100,-0.1,x[l-1],1.0);
        // data
        graph_serieline(graph, data_xy(x, y[0], 0, x[l - 1], 0, 1.0), 3, 2, new int[]{255, 0, 0}, new int[]{128, 0, 0});
        graph_serieline(graph, data_xy(x, y[1], 0, x[l - 1], 0, 1.0), 3, 2, new int[]{255, 0, 0}, new int[]{0, 128, 0});
        graph_serieline(graph, data_xy(x, y[2], 0, x[l - 1], 0, 1.0), 3, 2, new int[]{255, 0, 0}, new int[]{0, 0, 128});
        graph_serieline(graph, data_xy(x, y[3], 0, x[l - 1], 0, 1.0), 3, 2, new int[]{255, 0, 0}, new int[]{128, 128, 0});
        graph = graph_axis(graph, 50, 50, 50, 50);
        colorarr_show(graph);
    }


    public static double[] likelihood_s(int s) {
        // new
        denominator = 0;
        for (int i = 0, il = numerator.length; i < il; i++) {
            numerator[i] = 0;
        }

        // initial sample
        System.out.println("===================== iteration of " + s + " times ! =====================");
        for (int i = 0; i < s; i++) {
            likelihood_a();
        }
        double[] p_likelihood = new double[numerator.length];
        for (int i = 0, il = numerator.length; i < il; i++) {
            p_likelihood[i] = numerator[i] / denominator;
            System.out.println((i + 1) + "___" + p_likelihood[i]);
        }
        return p_likelihood;
    }

    public static double denominator = 0;
    public static double[] numerator = new double[10];
    public static double[] p_estimate = new double[10];
    public static int[] whichB = new int[]{2, 5, 8, 10}; //
    public static int[] B = new int[10];
    public static double alpha = 0.1;
    public static double coeffienct = (1 - alpha) / (1 + alpha);
    public static int Z = 128;

    public static void likelihood_a() {
        // Bi~P(Bi=1)=1/2
        // fB=2*(i-1)
        int FB = 0;
        for (int i = 9; i >= 0; i--) {
            FB *= 2;
            if (Math.random() > 0.5)
                B[i] = 1;
            else
                B[i] = 0;
            //System.out.print(B[i]);
            FB += B[i];
        }
        //System.out.println("____"+FB);


        double PQE = coeffienct * Math.pow(alpha, Math.abs(Z - FB));

        denominator += PQE;
        for (int i = 0, il = whichB.length; i < il; i++) {
            int index = whichB[i] - 1;
            if (B[index] == 1) {
                numerator[index] += PQE;
            }
        }
    }


    public static int[][] data_xy(int[] x, double[] y) {
        int l = x.length;
        int max_x = x[l - 1];
        int min_x = x[0];
        double max_y = y[0];
        double min_y = y[0];
        for (int i = 0; i < l; i++) {
            if (max_y < y[i])
                max_y = y[i];
            if (min_y > y[i])
                min_y = y[i];
        }

        return data_xy(x, y, (double) min_x, (double) max_x, min_y, max_y);
    }

    public static int[][] data_xy(int[] x, double[] y, double x1, double x2, double y1, double y2) {
        int l = x.length;
        double datx = x2 - x1;
        double daty = y2 - y1;
        int[][] xy = new int[2][l];
        for (int i = 0; i < l; i++) {
            xy[0][i] = (int) ((x[i] - x1) * g_w / datx);
            xy[1][i] = (int) ((y[i] - y1) * g_h / daty);
        }
        return xy;
    }


    public static void graph_serieline(int[][][] color, int[][] xy) {
        graph_serieline(color, xy, 3, 1, new int[]{255, 0, 0}, new int[]{0, 255, 0});
    }

    public static void graph_serieline(int[][][] color, int[][] xy, int r_point, int r_line, int[] color_point, int[] color_line) {
        int[] x = xy[0];
        int[] y = xy[1];
        int l = x.length;
        // line
        for (int i = 0; i + 1 < l; i++) {
            //System.out.println("==============================");
            //System.out.println(x[i]+"___"+y[i]+"__________"+x[i+1]+"___"+y[i+1]);
            graph_line(color, x[i], y[i], x[i + 1], y[i + 1], r_line, color_line);
        }
        // point
        for (int i = 0; i < l; i++) {
            graph_pointxy(color, x[i], y[i], r_point, color_point);
        }
    }

    public static int[] double_xyint(double x, double y, double x1, double y1, double x2, double y2) {
        System.out.println(x1 + "__" + y1 + "__" + x2 + "__" + y2 + "_________" + x + "___" + y);
        System.out.println((int) ((x - x1) * g_w / (x2 - x1)) + "___" + (int) ((y - y1) * g_h / (y2 - y1)));
        return new int[]{(int) ((x - x1) * g_w / (x2 - x1)), (int) ((y - y1) * g_h / (y2 - y1))};
    }

    public static void graph_line_relative(int[][][] color, double x1, double y1, double x2, double y2, int r, int[] linecolor, double xt1, double yt1, double xt2, double yt2) {
        int[] point1 = double_xyint(x1, y1, xt1, yt1, xt2, yt2);
        int[] point2 = double_xyint(x2, y2, xt1, yt1, xt2, yt2);
        System.out.println(point1[0] + "___" + point1[1] + "___" + point2[0] + "___" + point2[1]);
        graph_line(color, point1[0], point1[1], point2[0], point2[1], r, linecolor);
    }

    public static void graph_line(int[][][] color, int x1, int y1, int x2, int y2, int r, int[] linecolor) {
        int datx = x2 - x1;
        int daty = y2 - y1;
        if (datx == 0 && daty == 0) { // single point
            graph_pointxy(color, x1, y1, r, linecolor);
        } else if (Math.abs(datx) >= Math.abs(daty)) { // according to x
            if (datx > 0) {
                for (int i = x1; i <= x2; i++) {
                    //System.out.println(i+"___"+(y1+i*(daty)/datx));
                    graph_pointxy(color, i, y1 + (i - x1) * (daty) / datx, r, linecolor);
                }
            } else
                graph_line(color, x2, y2, x1, y1, r, linecolor);
        } else { // according to y
            if (daty > 0) {
                for (int j = y1; j <= y2; j++) {
                    graph_pointxy(color, x1 + (j - y1) * datx / daty, j, r, linecolor);
                }
            } else
                graph_line(color, x2, y2, x1, y1, r, linecolor);
        }
    }

    public static void graph_pointxy(int[][][] color, int x, int y, int r, int[] pointcolor) {
        int w = color.length;
        int h = color[0].length;
        for (int i = x - r; i <= x + r; i++) {
            if (i >= 0 && i < w)
                for (int j = y - r; j <= y + r; j++) {
                    if (j >= 0 && j < h)
                        for (int k = 0; k < 3; k++) {
                            color[i][j][k] = pointcolor[k];
                        }
                }
        }
    }


    public static int g_w = 0;
    public static int g_h = 0;

    public static int[][][] graph_new(int w, int h, int[] backgd) {
        g_w = w;
        g_h = h;
        int[][][] color = new int[g_w][g_h][3];
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                for (int k = 0; k < 3; k++) {
                    color[i][j][k] = backgd[k];
                }
            }
        }
        return color;
    }


    // last step !!!
    public static int[][][] graph_axis(int[][][] color, int left, int right, int down, int up) {
        int w = color.length;
        int h = color[0].length;
        int[][][] color_new = new int[w + left + right][h + up + down][];
        // graph
        for (int i = left, il = w + left; i < il; i++) {
            for (int j = up, jl = up + h; j < jl; j++) {
                //System.out.println((i-left)+"___"+(h-(j-up)));
                color_new[i][j] = color[i - left][h - 1 - (j - up)];
            }
        }

        // other
        // left
        int[] backgd = new int[]{100, 100, 100};
        for (int i = 0; i < left; i++) {
            for (int j = 0, jl = h + up + down; j < jl; j++) {
                color_new[i][j] = backgd;
            }
        }
        // right
        for (int i = w + left, il = w + left + right; i < il; i++) {
            for (int j = 0, jl = h + up + down; j < jl; j++) {
                color_new[i][j] = backgd;
            }
        }
        // up and down
        for (int i = left, il = left + w; i < il; i++) {
            for (int j = 0; j < up; j++) {
                color_new[i][j] = backgd;
            }
            for (int j = up + h, jl = up + h + down; j < jl; j++) {
                color_new[i][j] = backgd;
            }
        }
        return color_new;
    }

    public static void colorarr_show(int[][][] color) throws IOException {
        pixelarr_show(colorarr_pixelarr(color));
    }

    public static void pixelarr_show(int[][] pixel) throws IOException {
        UUID uuid = UUID.randomUUID();
        String plotpath = uuid.toString() + ".jpg"; // imagedir + "\\" +
        pixelarr_path(pixel, plotpath);
        path_open(plotpath);
    }

    public static void path_open(String path) throws IOException {
        Runtime.getRuntime().exec("rundll32 url.dll FileProtocolHandler " + path);
    }

    public static void pixelarr_path(int[][] pixel, String path) {
        int lx = pixel.length;
        int ly = pixel[0].length;
        BufferedImage img = new BufferedImage(lx, ly, BufferedImage.TYPE_INT_BGR);

        for (int i = 0; i < lx; i++)
            for (int j = 0; j < ly; j++) {
                img.setRGB(i, j, pixel[i][j]);
            }

        try {
            ImageIO.write(img, "bmp", new File(path));// jpeg may lose some information; bmp
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static int[][] colorarr_pixelarr(int[][][] color) {
        int lx = color.length;
        int ly = color[0].length;
        int[][] pixel = new int[lx][ly];
        for (int i = 0; i < lx; i++) {
            for (int j = 0; j < ly; j++) {
                //System.out.println(i+"__"+j+"__"+color[i][j][0]+"__"+color[i][j][1]+"__"+color[i][j][2]);
                Color c = new Color(color[i][j][0], color[i][j][1], color[i][j][2]);
                pixel[i][j] = c.getRGB();
            }
        }

        return pixel;
    }
}

