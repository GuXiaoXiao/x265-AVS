/*****************************************************************************
 * Copyright (C) 2013 x265 project
 *
 * Authors: Steve Borho <steve@borho.org>
 *          Min Chen <chenm003@163.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02111, USA.
 *
 * This program is also available under a commercial proprietary license.
 * For more information, contact us at license @ x265.com.
 *****************************************************************************/

#ifndef X265_MOTIONESTIMATE_H
#define X265_MOTIONESTIMATE_H

#include "primitives.h"
#include "reference.h"
#include "mv.h"
#include "bitcost.h"
#include "yuv.h"

namespace X265_NS {
// private x265 namespace

class MotionEstimate : public BitCost
{
protected:

    intptr_t blockOffset;					// 搜索块首地址相对于帧首地址的偏移量  
    
    int ctuAddr;
    int absPartIdx;  // part index of PU, including CU offset within CTU

    int searchMethod;					// ME搜索算法  
    int subpelRefine;					// subme强度  

	int blockwidth;						// 当前搜索块的宽度
	int blockheight;					// 在x265框架中暂时没有任何用途  

    pixelcmp_t sad;						// 用于计算块的SAD值 (src, srcStride, dst, dstStride,);  
    pixelcmp_x3_t sad_x3;				// 同时计算3个MV对应的3个SAD值			template<int lx, int ly> void sad_x3(const pixel* pix1, const pixel* pix2, const pixel* pix3, const pixel* pix4, intptr_t frefstride, int32_t* res)  
    pixelcmp_x4_t sad_x4;				// 同时计算4个MV对应的4个SAD值 template<int lx, int ly> oid sad_x4(const pixel* pix1, const pixel* pix2, const pixel* pix3, const pixel* pix4, const pixel* pix5, intptr_t frefstride, int32_t* res)  
    pixelcmp_ads_t ads;					// 计算SATD值，计算过程可以查看satd_4x4函数(fencIntra, cuSize, prediction, cuSize)
    pixelcmp_t satd;
    pixelcmp_t chromaSatd;

    MotionEstimate& operator =(const MotionEstimate&);

public:

    static const int COST_MAX = 1 << 28;

    uint32_t* integral[INTEGRAL_PLANE_NUM];
    Yuv fencPUYuv;						// 待搜索块的缓存,大小为64x64，将来搜索块会先copy到此缓存  
    int partEnum;
    bool bChromaSATD;					// 是否计算chroma分量的satd。只有在subpelRefine大于2时，在分像素ME时才会计算chroma的satd  

	/** 函数功能             ：初始化ME，searchMethod 默认hex，subme 默认2
	* \返回                  ：null * */
    MotionEstimate();
    ~MotionEstimate();

    static void initScales();
    static int hpelIterationCount(int subme);

	/** 函数功能   ： 初始化搜索算法、创建待搜索块的缓存
	/*\参数  method： 搜索方法
	/*\参数  refine： subme强度
	* \参数     csp:  图像格式
	* \返回        ： null */
    void init(int csp);

    /* Methods called at slice setup */
	/* Methods called at slice setup */
	/** 函数功能   ： 设置me对应的asm函数，copy待搜索块数据到待搜索块的缓存
	/*\参数   fencY： 当前编码帧的帧首地址
	/*\参数  stride： 原始帧步长
	/*\参数  offset： 当前搜索块首地址相对于帧首地址的偏移量
	/*\参数  pwidth： 当前搜索块的宽度
	/*\参数 pheight： 当前搜索块的高度
	* \返回        ： null */
    void setSourcePU(pixel *fencY, intptr_t stride, intptr_t offset, int pwidth, int pheight, const int searchMethod, const int subpelRefine);
    void setSourcePU(const Yuv& srcFencYuv, int ctuAddr, int cuPartIdx, int puPartIdx, int pwidth, int pheight, const int searchMethod, const int subpelRefine, bool bChroma);

    /* buf*() and motionEstimate() methods all use cached fenc pixels and thus
     * require setSourcePU() to be called prior. */
	/** 函数功能   ： 计算当前块与参考块的SATD值
	/*\参数    fref： 参考块首地址
	/*\参数  stride： 参考块步长步长
	* \返回        ： 当前块与参考块的SATD值 */
    inline int bufSAD(const pixel* fref, intptr_t stride)  { return sad(fencPUYuv.m_buf[0], FENC_STRIDE, fref, stride); }

    inline int bufSATD(const pixel* fref, intptr_t stride) { return satd(fencPUYuv.m_buf[0], FENC_STRIDE, fref, stride); }

    inline int bufChromaSATD(const Yuv& refYuv, int puPartIdx)
    {
        return chromaSatd(refYuv.getCbAddr(puPartIdx), refYuv.m_csize, fencPUYuv.m_buf[1], fencPUYuv.m_csize) +
               chromaSatd(refYuv.getCrAddr(puPartIdx), refYuv.m_csize, fencPUYuv.m_buf[2], fencPUYuv.m_csize);
    }

    int motionEstimate(ReferencePlanes* ref, const MV & mvmin, const MV & mvmax, const MV & qmvp, int numCandidates, const MV * mvc, int merange, MV & outQMv, pixel *srcReferencePlane = 0);

    int subpelCompare(ReferencePlanes* ref, const MV &qmv, pixelcmp_t);

protected:
	/**
	* 函数功能：星形ME搜索
	/*  调用范围        只在MotionEstimate::motionEstimate函数中被调用
	* \参数  ref        参考帧
	* \参数  mvmin      输出的实际搜索范围（左边界和上边界）
	* \参数  mvmax      输出的实际搜索范围（下边界和右边界）
	* \参数  bmv        从AMVP得到的预测MV，并返回最优的MV
	* \参数  bcost      预测MV对应的cost，并返回最优的cost
	* \参数  bPointNr   返回最优的MV对应的位置标号，该位置标号在下面ME的搜索模板中标出
	* \参数  bDistance  返回最优的MV对应的步长
	* \参数  earlyExitIters 提前跳出的迭代次数
	* \参数  merange    输入的ME搜索范围
	*/
    inline void StarPatternSearch(ReferencePlanes *ref,
                                  const MV &       mvmin,
                                  const MV &       mvmax,
                                  MV &             bmv,
                                  int &            bcost,
                                  int &            bPointNr,
                                  int &            bDistance,
                                  int              earlyExitIters,
                                  int              merange);
};
}

#endif // ifndef X265_MOTIONESTIMATE_H
