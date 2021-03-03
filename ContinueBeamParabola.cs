using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Autodesk.AutoCAD.Geometry2;

namespace TX_TOOLBOX
{
    class ContinueBeamParabolaCurve : ICloneable
    {
        ////目前处理根据原GetDownCurves_Parabola函数处理，并不取出实际位置线型
        ////分段长度---从跨中开始向墩顶(原点是跨中的点)
        ////------------------求方程曲线(只有曲线其他不包括)-----------------
        ////纯抛物线曲线--起点是梁根部位置-向上向右
        ////变为梁底线还要考虑跨中平段
        ////       
        ////        Y                     3   2    1    0   
        ////       /|\               4    ________________
        ////        |            5 ______/               (跨中)
        ////        |       6  ___/           
        ////        |   ...___/
        ////        | n___/  
        ////  n+1   | /
        ////  ------|------------------------------------->X
        ////   (墩顶)   (0,0)     y = a * x^n
        ////                

        public ContinueBeamParabolaCurve()
        {
            MaxFloorThick = 500;
            MinFloorThick = 250;
            LengthBottom = 20000;
            LengthFloor = 20000;
            LengthCrossBeam = 1000;
            LengthHoriSup = 1500;
            PowOfBottom = 2.0;
            PowOfFloor = 2.0;
            CvFloor = 500;
            MaxBeamHeight = 3000;
            MinBeamHeight = 2000;
        }
        public ContinueBeamParabolaCurve(SpanLK span, bool isSpanLeft)
        {
            IsByVirtualPoint = span.IsByVPoint;
            MaxFloorThick = isSpanLeft ? span.m_abt1 : span.m_abt2;
            MinFloorThick = span.m_bt;
            LengthBottom = isSpanLeft ? span.DownCurveL1 : span.DownCurveL2;
            LengthFloor = isSpanLeft ? span.m_DL_chge1 : span.m_DL_chge2;
            LengthCrossBeam = isSpanLeft ? span.m_HL1 : span.m_HL2;
            LengthHoriSup = isSpanLeft ? span.DownCurveX1 : span.DownCurveX2;
            PowOfBottom = span.DownCurvePow;
            PowOfFloor = span.DibanChangePow;
            CvFloor = isSpanLeft ? span.cvXd1 : span.cvXd2;
            MaxBeamHeight = isSpanLeft ? span.Hmax1 : span.Hmax2;
            MinBeamHeight = span.Hm;
        }

        private bool isByVirtualPoint = false;
        /// <summary>
        /// 底板加厚厚度
        /// </summary>
        public double MaxFloorThick { get; set; }
        /// <summary>
        /// 底板标准厚度
        /// </summary>
        public double MinFloorThick { get; set; }
        /// <summary>
        /// 底缘抛物线水平长度
        /// </summary>
        public double LengthBottom { get; set; }
        /// <summary>
        /// 底板加厚段长度
        /// </summary>
        public double LengthFloor { get; set; }
        /// <summary>
        /// 横梁宽度
        /// </summary>
        public double LengthCrossBeam { get; set; }
        /// <summary>
        /// 底缘水平段长度
        /// </summary>
        public double LengthHoriSup { get; set; }
        /// <summary>
        /// 底缘抛物线次幂
        /// </summary>
        public double PowOfBottom { get; set; }
        /// <summary>
        /// 底板加厚次幂
        /// </summary>
        public double PowOfFloor { get; set; }
        /// <summary>
        /// 横梁处加厚倒角长度
        /// </summary>
        public double CvFloor { get; set; }
        /// <summary>
        /// 墩顶梁高
        /// </summary>
        public double MaxBeamHeight { get; set; }
        /// <summary>
        /// 跨中梁高
        /// </summary>
        public double MinBeamHeight { get; set; }
        /// <summary>
        /// 虚交点
        /// </summary>
        public bool IsByVirtualPoint
        {
            get { return isByVirtualPoint && LengthCrossBeam == LengthHoriSup && LengthFloor == LengthBottom; }
            set { isByVirtualPoint = value; }
        }

        public ParabolaEquation BottomEquation
        {
            get
            {
                Point3d apex = new Point3d(LengthBottom, MaxBeamHeight - MinBeamHeight, 0);
                ParabolaEquation equation = ParabolaEquation.CreateEquation(PowOfBottom, apex, new Point3d());
                return equation;
            }
        }

        public object Clone()
        {
            return this.MemberwiseClone();
        }

        public double GetBottomStretchedPow(double factorX, double factorY)
        {
            double pow = ParabolaEquation.GetStretchedPow(this.BottomEquation, factorX, factorY);
            return pow;
        }

        public double GetFloorStretchedPow(double factorX, double factorY)
        {
            if (IsByVirtualPoint)
            {
                Point3d apex = new Point3d(LengthBottom, MaxBeamHeight - MinBeamHeight + MinFloorThick, 0);
                Point3d pt2 = new Point3d(0, MaxFloorThick, 0);
                ParabolaEquation equationFloor = ParabolaEquation.CreateEquation(PowOfFloor, apex, pt2);
                return ParabolaEquation.GetStretchedPow(equationFloor, factorX, factorY);
            }
            else
            {
                Point3d apexFloor = new Point3d(LengthFloor, 0, 0);
                Point3d pt2Floor = new Point3d(0, MaxFloorThick - MinFloorThick, 0);
                ParabolaEquation equationFloor = ParabolaEquation.CreateEquation(PowOfFloor, apexFloor, pt2Floor);
                return ParabolaEquation.GetStretchedPow(equationFloor, factorX, factorY);
            }
        }

        public ContinueBeamParabolaCurve GetStretchedParabolaCurve(double factorX, double factorY)
        {
            ContinueBeamParabolaCurve newParabolaCurve = this.Clone() as ContinueBeamParabolaCurve;
            newParabolaCurve.LengthBottom *= factorX;
            newParabolaCurve.LengthFloor *= factorX;
            newParabolaCurve.MaxBeamHeight *= factorY;
            newParabolaCurve.MinBeamHeight *= factorY;
            return newParabolaCurve;
        }

        /// <summary>
        /// 这个函数线的方向和现有的还不太一样，还不能直接替换  TODO:修改点位顺序
        /// </summary>
        public TxPolyline GetBottomPolyline(double precision)
        {
            double lengthFloorNormal = Math.Max(0, LengthHoriSup + LengthBottom - LengthCrossBeam - LengthFloor);                              //变高范围内底板等厚段长度,有可能是会小于0，但是目前没有去支持这样的情况，直接判定掉了
            double lengthFloorChange = LengthHoriSup + LengthBottom - lengthFloorNormal - Math.Max(LengthCrossBeam + CvFloor, LengthHoriSup);  //底板厚度变化段 
            double lengthFloorRemain = Math.Max(0, LengthBottom - lengthFloorNormal - lengthFloorChange);                                      //剩余倒角段部分，未伸至全长

            List<double> segFloorNormal = TX_Math.Design.BJ(lengthFloorNormal, precision, TX_Math.Design.BJType.Average);
            List<double> segFloorChange = TX_Math.Design.BJ(lengthFloorChange, precision, TX_Math.Design.BJType.Average);
            List<double> segFloorRemain = TX_Math.Design.BJ(lengthFloorRemain, precision, TX_Math.Design.BJType.Average);

            List<double> segCtl = new List<double>();
            segCtl.AddRange(segFloorRemain);
            segCtl.AddRange(segFloorChange);
            segCtl.AddRange(segFloorNormal);

            TxPolyline result = new TxPolyline();
            //添加墩顶水平段
            if (LengthHoriSup > 0)
                result.AddVertexAt(new Point3d(-LengthHoriSup, 0, 0));

            ParabolaEquation bottomEquation = this.BottomEquation;

            //从左往右顺序添加节点
            double tempX = 0;
            for (int i = 0; i < segCtl.Count; i++)
            {
                Point3d tempPt = bottomEquation.GetPoint(tempX);
                result.AddVertexAt(tempPt);
                tempX += segCtl[i];
            }

            //添加顶点
            result.AddVertexAt(BottomEquation.Apex);

            return result;
        }

        /// <summary>
        /// 这个函数线的方向和现有的还不太一样，还不能直接替换  TODO:修改点位顺序
        /// </summary>
        public TxPolyline GetFloorPolyline(double precision)
        {
            //一种是虚交点，这个应该是直接用方程计算
            //一种是算控制点分段计算
            if (IsByVirtualPoint)
            {
                double lengthFloorChange = LengthFloor - CvFloor;  //底板厚度变化段 
                List<double> segFloor = TX_Math.Design.BJ(lengthFloorChange, precision, TX_Math.Design.BJType.Average);

                Point3d apex = new Point3d(LengthBottom, MaxBeamHeight - MinBeamHeight + MinFloorThick, 0);
                Point3d pt2 = new Point3d(0, MaxFloorThick, 0);
                ParabolaEquation equation = ParabolaEquation.CreateEquation(PowOfFloor, apex, pt2);
                TxPolyline result = new TxPolyline();

                //从左往右顺序添加节点
                double tempX = CvFloor;
                for (int i = 0; i < segFloor.Count; i++)
                {
                    Point3d tempPt = equation.GetPoint(tempX);
                    result.AddVertexAt(tempPt);
                    tempX += segFloor[i];
                }
                //添加顶点
                result.AddVertexAt(apex);

                return result;
            }
            else
            {
                double lengthFloorNormal = Math.Max(0, LengthHoriSup + LengthBottom - LengthCrossBeam - LengthFloor);                        //变高范围内底板等厚段长度,有可能是会小于0，但是目前没有去支持这样的情况，直接判定掉了
                double lengthFloorChange = LengthHoriSup + LengthBottom - lengthFloorNormal - Math.Max(LengthCrossBeam + CvFloor, LengthHoriSup);  // 底板厚度变化段 

                List<double> segFloorNormal = TX_Math.Design.BJ(lengthFloorNormal, precision, TX_Math.Design.BJType.Average);
                List<double> segFloorChange = TX_Math.Design.BJ(lengthFloorNormal, precision, TX_Math.Design.BJType.Average);

                List<double> segCtl = new List<double>();
                segCtl.AddRange(segFloorChange);
                segCtl.AddRange(segFloorNormal);

                //底板上缘看原先的书写是底板厚度的变化值是按照抛物线来计算的，与底板下缘线型采用抛物线有一些不一样
                Point3d apexFloor = new Point3d(LengthFloor, 0, 0);
                Point3d pt2Floor = new Point3d(0, MaxFloorThick - MinFloorThick, 0);
                ParabolaEquation equationFloor = ParabolaEquation.CreateEquation(PowOfFloor, apexFloor, pt2Floor);
                TxPolyline result = new TxPolyline();

                //这里有个异常情况是倒角段也没伸出平段，这里特殊处理一下添加一个平段
                if (LengthCrossBeam + CvFloor - LengthHoriSup < 0)
                    result.AddVertexAt(new Point3d(LengthCrossBeam + CvFloor - LengthHoriSup, MaxFloorThick, 0));

                //从左往右顺序添加节点
                double tempX = Math.Max(0, LengthCrossBeam + CvFloor - LengthHoriSup);
                for (int i = 0; i < segCtl.Count; i++)
                {
                    Point3d tempPt = BottomEquation.GetPoint(tempX);
                    double thick = equationFloor.GetY(tempX);
                    tempPt += Vector3d.YAxis * thick;
                    result.AddVertexAt(tempPt);
                    tempX += segCtl[i];
                }

                //添加顶点
                result.AddVertexAt(BottomEquation.Apex + Vector3d.YAxis * MinFloorThick);

                return result;
            }
        }

    }

    struct ParabolaEquation
    {
        private ParabolaEquation(double a, double h, double k, double pow)
        {
            this.a = a;
            this.h = h;
            this.k = k;
            this.pow = pow;
            this.factorX = 1.0;
            this.factorY = 1.0;
        }

        //现阶段只能用顶点式
        private double a;
        private double h;
        private double k;
        private double pow;
        private double factorX;
        private double factorY;

        public double Pow { get { return pow; } }
        public Point3d Apex { get { return new Point3d(h, k, 0); } }

        public double GetY(double x)
        {
            double calX = x / factorX;
            double y = a * Math.Pow(Math.Abs(calX - h), pow) + k;
            return factorY * y;
        }

        public Point3d GetPoint(double x)
        {
            Point3d pt = new Point3d(x, GetY(x), 0);
            return pt;
        }

        /// <summary>
        /// 提供给现有的调用的，主要这部分都写好了，单独提取感觉有点麻烦
        /// </summary>
        /// <param name="equation"></param>
        /// <param name="factorX"></param>
        /// <param name="factorY"></param>
        /// <returns></returns>
        public static double GetStretchedPow(ParabolaEquation equation, double factorX, double factorY)
        {
            ParabolaEquation newEquation = ParabolaEquation.StretchParabola(equation, factorX, factorY);
            return newEquation.Pow;
        }

        public static ParabolaEquation StretchParabola(ParabolaEquation equation, double factorX, double factorY)
        {
            ParabolaEquation newEquation = equation;
            //关键在对于次幂的修正，这里做的是以原点为中心的拉伸，因为是顶点式，所以只需要修改顶点就能做到平移
            //要计算次幂的修正值，需要除顶点外的另外两个不对称点，这边就直接从方程中随便取两个点了
            Point3d apex = equation.Apex;
            Point3d pt1 = equation.GetPoint(equation.h + 1000);
            Point3d pt2 = equation.GetPoint(equation.h + 2000);
            //点坐标计算上因数，即为拉伸后的曲线上坐标点
            apex = new Point3d(apex.X * factorX, apex.Y * factorY, 0);
            pt1 = new Point3d(pt1.X * factorX, pt1.Y * factorY, 0);
            pt2 = new Point3d(pt2.X * factorX, pt2.Y * factorY, 0);

            //方程转换知 pow = log(y-k)-log(a) 底数为(x-h)
            //则前半部分的对数取值log(y-k)
            double lp1 = Math.Log(Math.Abs(pt1.Y - apex.Y), pt1.X - apex.X);
            double lp2 = Math.Log(Math.Abs(pt2.Y - apex.Y), pt2.X - apex.X);
            double lnX1 = Math.Log(Math.Abs(pt1.X - apex.X));
            double lnX2 = Math.Log(Math.Abs(pt2.X - apex.X));
            double lna = lnX1 * lnX2 * (lp1 - lp2) / (lnX2 - lnX1);

            double a = Math.Sign(pt1.Y - apex.Y) * Math.Pow(Math.E, lna);
            double h = apex.X;
            double k = apex.Y;
            double pow = Math.Round(Math.Log((pt1.Y - k) / a, pt1.X - h), 2);      //纵向里面的修正应该就采用这个去修正就可以对上了
            return ParabolaEquation.CreateEquation(a, h, k, pow);
        }

        public static ParabolaEquation CreateEquation(double a, double h, double k, double pow)
        {
            return new ParabolaEquation(a, h, k, pow);
        }

        public static ParabolaEquation CreateEquation(double pow, Point3d apex, Point3d pt2)
        {
            double h = apex.X;
            double k = apex.Y;
            double a = (pt2.Y - k) / Math.Pow(Math.Abs(pt2.X - h), pow);
            return new ParabolaEquation(a, h, k, pow);
        }

    }
}
