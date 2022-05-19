import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.time.LocalDateTime;
import java.time.ZoneId;
import java.util.ArrayList;
import java.util.Date;
import java.util.Scanner;

import org.ujmp.core.DenseMatrix;
import org.ujmp.core.Matrix;

import javax.transaction.xa.XAResource;


class satelitepos {
    double Mk;

    double Ek;
    private double dtr;//相对论改正数
    protected final double c=2.99792458E08;
    String address1;//txt文件位置
    String satellite;//PRN号
    protected double IDOT;
    protected double sqrt_a;
    protected double IODE;
    protected double delta_n;
    protected double t;
    protected double toe;
    protected double M0;
    protected double e;
    protected double omg;
    protected double Cuc;
    protected double Cus;
    protected double Crc;
    protected double Crs;
    protected double Cic;
    protected double Cis;
    protected double i0;
    protected double OMG0;
    protected double OMG0_DOT;
    protected double [][]position=new double[1][3];
    double dt;
    ArrayList<String> arr2=new ArrayList<String>();
    ArrayList<String> arr3=new ArrayList<String>();
    double deltat;
    double Rs;//是卫星的几何距离



    Matrix Xr=DenseMatrix.Factory.zeros(3,1);//测站点的概略坐标
    double [][]poschange=new double [1][3];
    Matrix Pc=DenseMatrix.Factory.importFromArray(poschange);//存储在经过地球自转改正后的卫星坐标矩阵



    public void calculateLocation() {
        //长半径
        double a = Math.pow(this.sqrt_a, 2);
        //地球重力常数
        double constGM = 3.986004418e14;
        //平均角速度
        double n0 = Math.sqrt(constGM / Math.pow(a, 3));
        double n = n0 + this.delta_n;
        //平近点角
        Double tk = this.t - this.toe;
        this.Mk = this.M0 + n * (this.t - this.toe);
        //偏近点角
        double Ek0 = 0;
        double Ek1 = Mk;
        while (Math.abs(Ek1 - Ek0) > 2e-12) {
            Ek0 = Ek1;
            Ek1 = Mk + this.e * Math.sin(Ek0);
        }
        this.Ek = Ek0;
        //真近点角
        double vk = Math.atan((Math.pow(1-e*e,0.5)*Math.sin(Ek))/(Math.cos(Ek)-this.e));
        //升交角距
        double uk = vk + this.omg;
        //扰动改正
        double delta_uk = this.Cuc * Math.cos(2 * uk) + this.Cus * Math.sin(2 * uk);
        double delta_rk = this.Crc * Math.cos(2 * uk) + this.Crs * Math.sin(2 * uk);
        double delta_ik = this.Cic * Math.cos(2 * uk) + this.Cis * Math.sin(2 * uk);
        //赤经 半径 轨道倾角改正
        double u = uk + delta_uk;
        double r = a * (1 - this.e * Math.cos(Ek)) + delta_rk;
        double i = this.i0 + delta_ik+this.IDOT*tk;
        //卫星在轨道上的位置
        double x = r * Math.cos(u);
        double y = r * Math.sin(u);
        //地球自转角速度
        double omge = 7.292115e-5;
        //历元t升交点的赤纬
        double lambda = this.OMG0 + this.OMG0_DOT* tk - omge * this.toe;
        //卫星在地心地固坐标系中的位置
        double x_E=x*Math.cos(lambda)-y*Math.cos(i)*Math.sin(lambda);
        double y_E=x*Math.sin(lambda)+y*Math.cos(i)*Math.cos(lambda);
        double z_E=y*Math.sin(i);
        //卫星在BDCS坐标系的坐标
        double fi=omge*tk;
        double five=180/Math.PI*5;
        double x_R=Math.cos(fi)*x_E+Math.sin(fi)*Math.cos(five)*y_E-Math.sin(fi)*Math.sin(five)*z_E;
        double y_R=-1*Math.sin(fi)*x_E+Math.cos(fi)*Math.cos(five)*y_E-Math.sin(five)*Math.cos(fi)*z_E;
        double z_R=Math.sin(five)*y_E+Math.cos(five)*z_E;
        double distance=Math.pow(x_R,2)+Math.pow(y_R,2)+Math.pow(z_R,2);
        //System.out.println("距离地面高度为："+(Math.sqrt(distance)-6378137));
        this.position[0][0]=x_E;
        this.position[0][1]=y_E;
        this.position[0][2]=z_E;


    }
    public void getPOS (double []pos){
        for(int i=0;i<3;i++){
            pos[i]=this.position[0][i];
        }
    }



    public void read(){
        ArrayList<String> arr=new ArrayList<String>();
        ArrayList<String> arr1=new ArrayList<String>();
        String str[]=new String[36];
        try{

            Scanner sc=new Scanner(new File(address1));
            while (sc.hasNextLine()) {
                str = sc.nextLine().split(" ");

                for(int i=0;i<str.length;i++){

                    arr.add(str[i].replaceAll("D","E"));
                }
            }
        }catch(FileNotFoundException e){
            e.printStackTrace();
        }
        for (int i = 0; i <arr.size(); i++) {
            if (arr.get(i)!=null&&!arr.get(i).equals("")){
                arr1.add(arr.get(i));
            }
        }
        this.satellite=arr1.get(0);
        this.IODE=Double.parseDouble(arr1.get(10));
        //System.out.println("这是IODE："+this.IODE);
        //Crs
        this.Crs=Double.parseDouble(arr1.get(11));
        //System.out.println("这是Crs："+this.Crs);
        //角速度改正数
        this.delta_n=Double.parseDouble(arr1.get(12));
        //参考时刻的平近点角
        this.M0=Double.parseDouble(arr1.get(13));
        //Cuc
        this.Cuc=Double.parseDouble(arr1.get(14));
        //偏心率
        this.e=Double.parseDouble(arr1.get(15));
        //Cus
        this.Cus=Double.parseDouble(arr1.get(16));
        //长半轴的开方
        this.sqrt_a=Double.parseDouble(arr1.get(17));
        //toe
        this.toe=Double.parseDouble(arr1.get(18));
        //Cic
        this.Cic=Double.parseDouble(arr1.get(19));
        //参考时刻的升交点赤经
        this.OMG0=Double.parseDouble(arr1.get(20));
        //Cis
        this.Cis=Double.parseDouble(arr1.get(21));
        //参考时刻的轨道倾角
        this.i0=Double.parseDouble(arr1.get(22));
        //Crc
        this.Crc=Double.parseDouble(arr1.get(23));
        //参考时刻的近地点角距
        this.omg=Double.parseDouble(arr1.get(24));
        //升交点赤经变化率
        this.OMG0_DOT=Double.parseDouble(arr1.get(25));
        //轨道倾角变化率
        this.IDOT=Double.parseDouble(arr1.get(26));

        //toc 时钟时间的周内秒
        Integer year=Integer.parseInt(arr1.get(1));
        Integer month=Integer.parseInt(arr1.get(2));
        Integer date=Integer.parseInt(arr1.get(3));
        Integer hour=Integer.parseInt(arr1.get(4));
        Integer minu=Integer.parseInt(arr1.get(5));
        Integer sec=Integer.parseInt(arr1.get(6));
        LocalDateTime time=LocalDateTime.of(year, month, date, hour, minu, sec);
        Date dtTime=Date.from(time.atZone(ZoneId.systemDefault()).toInstant());
        Long timeSec=dtTime.getTime();

        LocalDateTime time0=LocalDateTime.of(1980,1,6,0,0,0);
        Date dtTime0=Date.from(time0.atZone(ZoneId.systemDefault()).toInstant());
        Long time0Sec=dtTime0.getTime();
        //求toc 时钟时间
        Long toc=((timeSec-time0Sec)/1000)%(604800);

        //钟差参数 a0,a1,a2
        Double a0=Double.parseDouble(arr1.get(7));
        Double a1=Double.parseDouble(arr1.get(8));
        Double a2=Double.parseDouble(arr1.get(9));
        //设置观测时间为每分钟的第13秒
        Long t0=toc;
        LocalDateTime time_ob=LocalDateTime.of(year, month, date, hour, minu, sec);
        this.dt=a0+a1*(t0-toc)+a2*Math.pow((t0-toc), 2);
        this.t=t0-dt;

        calculateLocation();



    }
    public double dis(double x1,double y1,double z1,double x2,double y2,double z2){
        double result=Math.pow(x1-x2,2)+Math.pow(y1-y2,2)+Math.pow(z1-z2,2);
        return Math.pow(result,0.5);
    }
    public void read(String address2){
        String str[]=new String[36];
        try{

            Scanner sc=new Scanner(new File(address2));
            while (sc.hasNextLine()) {
                str = sc.nextLine().split(" ");

                for(int i=0;i<str.length;i++){

                    this.arr2.add(str[i].replaceAll("D","E"));
                }
            }
        }catch(FileNotFoundException e){
            e.printStackTrace();
        }
        for (int i = 0; i <arr2.size(); i++) {
            if (this.arr2.get(i)!=null&&!this.arr2.get(i).equals("")){
                this.arr3.add(arr2.get(i));
            }
        }

    }
    public void diedai1(int n){

        shijiangaizheng();//计算相对论改正数
        double tao0=(Double.parseDouble(arr3.get(n))/this.c)+this.dt;//确定起始tao0
        this.deltat=this.t-tao0;//确定起始Tsi
        //read();//确定起始卫星位置
        poschange=position;//确定起始位置
        this.t=this.deltat+this.dtr;
        calculateLocation();
        this.Pc=DenseMatrix.Factory.importFromArray(poschange);//导入矩阵进行初始化
        this.Pc= satpos(poschange);//计算地球自转改正

        double Rs=dis(Xr.getAsDouble(0,0),Xr.getAsDouble(1,0),Xr.getAsDouble(2,0),Pc.getAsDouble(0,0),Pc.getAsDouble(1,0),Pc.getAsDouble(2,0));
        double tao1=Rs/this.c;
        int cnt=0;


        while(Math.abs(tao0-tao1)>1E-7) {

            tao0=tao1;
            this.deltat-=tao0;//计算下一个时间

            this.t=this.deltat+this.dtr;
            calculateLocation();

            this.Pc= satpos(poschange);//计算地球自转改正
            this.Rs=dis(Xr.getAsDouble(0,0),Xr.getAsDouble(1,0),Xr.getAsDouble(2,0),Pc.getAsDouble(0,0),Pc.getAsDouble(1,0),Pc.getAsDouble(2,0));//计算几何距离
            tao1 = this.Rs / this.c;//计算新的tao1
            cnt++;



        }





    }
    public void shijiangaizheng(){
        this.dtr=2290*(this.Ek-this.Mk)*1E-9;
    }

    public Matrix satpos(double [][]position)  {
        //n为在列表中的位数（由于读取数据时没有使用数组进行位数的增加而导致的叠加数）
        //地球自转改正卫星运行时的坐标
        Matrix a =DenseMatrix.Factory.importFromArray(position);
        final double omega = 7.2921151467e-5;
        double[][] result = new double[1][3];
        Matrix re = DenseMatrix.Factory.importFromArray(result);
        double[][] func = {{Math.cos(omega * deltat), Math.sin(omega * deltat), 0}, {-Math.sin(omega * deltat), Math.cos(omega * deltat), 0}, {0, 0, 1}};//地球自传改正系数矩阵
        Matrix function = DenseMatrix.Factory.importFromArray(func);
        re = function.mtimes(a.transpose());



        return re;

    }
    public  double [] getcos() {
        //计算卫星Si方向余弦
        double b0s1 = (Xr.getAsDouble(0,0)- position[0][0]) / this.Rs;
        double b1s1 = (Xr.getAsDouble(1,0)- position[0][1]) / this.Rs;
        double b2s1 = (Xr.getAsDouble(2,0)- position[0][2]) / this.Rs;
        double b3s1 = 1.0;
        double result[]={b0s1,b1s1,b2s1,b3s1};
        return result;
    }
    //public double calDtrop(){
        //计算对流层延迟改正量

    //}
   // public double calDiono(){
     //   //计算电离层延迟改正量

    //}
    public double yushu(int n){
        double Dtrop=0;
        double Diono=0;
        return Double.parseDouble(arr3.get(4+2*n))-dis(Xr.getAsDouble(0,0),Xr.getAsDouble(1,0),Xr.getAsDouble(2,0),poschange[0][0],poschange[0][1],poschange[0][2])+this.c*this.dt-Diono-Dtrop;
    }
    public void setXr(Matrix x0){
        this.Xr=x0;
    }
    public satelitepos(String address){
        this.address1=address;
    }

}

public class position {

    public static void main(String[] args)throws FileNotFoundException,IOException {
        int cnt=0;
        double sat1[] = new double[3];
        double sat2[] = new double[3];
        double sat3[] = new double[3];
        double sat4[] = new double[3];
        Matrix X = DenseMatrix.Factory.zeros(4, 1);//原始计算值

        satelitepos l1 = new satelitepos("C:\\Users\\Administrator\\Desktop\\GPS卫星1.txt");
        l1.read();
        l1.getPOS(sat1);
        satelitepos l2 = new satelitepos("C:\\Users\\Administrator\\Desktop\\GPS卫星2.txt");
        l2.read();
        l2.getPOS(sat2);
        satelitepos l3 = new satelitepos("C:\\Users\\Administrator\\Desktop\\GPS卫星3.txt");
        l3.read();
        l3.getPOS(sat3);
        satelitepos l4 = new satelitepos("C:\\Users\\Administrator\\Desktop\\GPS卫星4.txt");
        l4.read();
        l4.getPOS(sat4);



        /*Matrix a2=DenseMatrix.Factory.importFromArray(sat1);
        Matrix b1=DenseMatrix.Factory.importFromArray(sat2);
        Matrix c=DenseMatrix.Factory.importFromArray(sat3);
        Matrix d=DenseMatrix.Factory.importFromArray(sat4);
        System.out.println(a2);
        System.out.println(b1);
        System.out.println(c);
        System.out.println(d);
        */

        l1.setXr(X);
        l2.setXr(X);
        l3.setXr(X);
        l4.setXr(X);
        l1.read("C:\\Users\\Administrator\\Desktop\\观测值提取.txt");
        l2.read("C:\\Users\\Administrator\\Desktop\\观测值提取.txt");
        l3.read("C:\\Users\\Administrator\\Desktop\\观测值提取.txt");
        l4.read("C:\\Users\\Administrator\\Desktop\\观测值提取.txt");

        Matrix result = DenseMatrix.Factory.zeros(4, 1);
        double[][] L = new double[1][4];



        l1.diedai1(4);//小迭代运行
        l2.diedai1(6);
        l3.diedai1(8);
        l4.diedai1(10);
        double yushu1 = l1.yushu(0);//四颗卫星计算余数项
        double yushu2 = l2.yushu(1);
        double yushu3 = l3.yushu(2);
        double yushu4 = l4.yushu(3);

        L[0][0] = yushu1;//存入数组
        L[0][1] = yushu2;
        L[0][2] = yushu3;
        L[0][3] = yushu4;

        double []A1 = l1.getcos();//计算方向余弦
        double []A2 = l2.getcos();
        double []A3 = l3.getcos();
        double []A4 = l4.getcos();
        double[][] a = {A1, A2, A3, A4};
        Matrix Azhen = DenseMatrix.Factory.importFromArray(a);
        Matrix Lzhen = DenseMatrix.Factory.importFromArray(L);
        Matrix a1 = Azhen.transpose().mtimes(Azhen);
        Matrix b = Azhen.transpose().mtimes(Lzhen.transpose());

        Matrix result1 = a1.mtimes(b);//改正数矩阵计算
        result = result.plus(result1);







        while (Math.abs(result.getAsDouble(0, 0) - X.getAsDouble(0, 0)) > 1E-3
                && Math.abs(result.getAsDouble(1, 0) - X.getAsDouble(1, 0)) > 1E-3
                && Math.abs(result.getAsDouble(2, 0) - X.getAsDouble(2, 0)) > 1E-3) {


            X = X.plus(result);//计算改正后坐标
            l1.setXr(X); //修改小迭代测站点坐标矩阵
            l2.setXr(X);
            l3.setXr(X);
            l4.setXr(X);
            l1.diedai1(4);//小迭代运行
            l2.diedai1(6);
            l3.diedai1(8);
            l4.diedai1(10);
            yushu1 = l1.yushu(0);//四颗卫星计算余数项
            yushu2 = l2.yushu(1);
            yushu3 = l3.yushu(2);
            yushu4 = l4.yushu(3);

            L[0][0] = yushu1;//存入数组
            L[0][1] = yushu2;
            L[0][2] = yushu3;
            L[0][3] = yushu4;

            A1 = l1.getcos();//计算方向余弦
            A2 = l2.getcos();
            A3 = l3.getcos();
            A4 = l4.getcos();
            double [][]a0 = {A1, A2, A3, A4};
            Azhen = DenseMatrix.Factory.importFromArray(a0);
            Lzhen = DenseMatrix.Factory.importFromArray(L);
            a1 = Azhen.transpose().mtimes(Azhen);
            b = Azhen.transpose().mtimes(Lzhen.transpose());

            result1 = a1.mtimes(b);//改正数矩阵计算
            result = result.plus(result1);
            cnt++;
            System.out.println(cnt);



        }


    }

}