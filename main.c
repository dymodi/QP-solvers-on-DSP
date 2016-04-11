//
 #include <stdio.h>
 #include <string.h>
 #include <stdlib.h>
 #include <math.h>
 #include <limits.h>
 #include <stdint.h>


 #include <xdc/std.h>
 #include <xdc/runtime/Error.h>
 #include <xdc/runtime/System.h>
 #include <xdc/runtime/Types.h>
 #include <xdc/runtime/Memory.h>
 #include <xdc/runtime/Timestamp.h>
 #include <ti/sysbios/BIOS.h>
 #include <ti/sysbios/knl/Task.h>

 #include "ti\platform\platform.h"
 #include "D:\Users\Work\Documents\MATLAB\ASM\WGS\wgsQP\wgsQP\mathlib_double.h"

 FILE *fp_Fx, *fp_cons_obey;


 //============================mpc_lib_DP.c==START======================================

 //============================mpc_lib_DP.c==END========================================


 //
 //Read seria port：
 int post_read_uart(uint8_t* msg, uint32_t msg_len,uint32_t delay)
 {
     uint32_t i;

	 for(i=0;i<msg_len;i++)
	 {
		 if(platform_uart_read(&(msg[i]),delay) != Platform_EOK)
		 {
			 return -1;
		 }
	 }
	 return 0;
 }

 //Write seria port：
 int post_write_uart( char* msg,uint32_t msg_len)
 {
     uint32_t i;

     for (i = 0; i < msg_len; i++)
     {
         if (platform_uart_write(msg[i]) != Platform_EOK)
         {
             return -1;
         }
     }
     return 0;
 }


 Void taskFxn(UArg a0, UArg a1)
 {
     printf("enter taskFxn()\n");
     int ndec,mc,wSize;
     int i, kk, solverOpt;
     int maxIter = 300;
     int maxmem;
     int nbc = 0; // All constraints are general constraints
	 int *w, *iterPoint;
     double Time_MPC_Total, Time_MPC_Worst, Time_MPC;
     double *H, *c, *A, *b, *H_ori, *c_ori, *invH, *invH_ori, *x_ini, *x, *xStar;
     double *Aw, *bw, *p, *gx, *lambda;
     double *tmpVec1, *tmpVec2, *tmpMat1, *tmpMat2, *tmpMat3, *tmpMat4, *tmpMat5, *tmpMat6;
     double *Av, *Af, *Lv, *Yv, *Rv, *pvStar, *uv, *vl, *wv, *P, *G, *lambdaf, *lambdal, *lx, *ux;
     int *orderPermu,*isInWl, *wf, *wl;	//Initial working set is set to be empty
     int hasCalloc = 0;

	 UInt32 start_time,end_time,dead_time;
	 UInt32 MPC_start_time,MPC_end_time,MPC_total_time;
	 Types_FreqHz freq;


	 //Communication Settings
	 char token_read = 255;
	 char token_write = 240;
	 char token_init = 230;
	 char token_succ = 250;
	 float tmpFloat1 = 0;
	 float tmpFloat2 = 0;
	 float iniFlag = -1;
	 float count_write = 1;
	 float count_read = 1;
	 float count_r_verify;
	 int t_f1, t_f2;
	 //uint32_t delay = 1700000;	//摸索出来的，不要随便乱改
	 uint32_t delay = 2000000;	//摸索出来的，不要随便乱改

	 uint8_t Buffer_read[12];
	 char Buffer_write[8];


	 Timestamp_getFreq(&freq);
	 //printf("freq = {0x%x : 0x%x}\n", freq.hi, freq.lo);
	 printf("Frequency is: %u\n",freq.lo);


	 for (kk = 0; kk < 1000; kk++) {
		 // 这里先读ndec, mc, transOpt
		 //读8字节数据
		 t_f1 = post_read_uart(&Buffer_read[0],16,delay);
		 memcpy(&count_r_verify,&Buffer_read[0],sizeof(float));
		 memcpy(&tmpFloat1,&Buffer_read[4],sizeof(float));
		 memcpy(&tmpFloat2,&Buffer_read[8],sizeof(float));
		 memcpy(&iniFlag,&Buffer_read[12],sizeof(float));

		 if(count_r_verify!=count_read) {
			 printf("Count_read Verify failed! \n\r");
			 printf("count_r_verify is %f:! \n\r",count_r_verify);
			 printf("count_read is %f:! \n\r",count_read);
			 printf("Is read successful? %d \n\r",t_f1);
			 printf("What's in the buffer? %c \n\r",Buffer_read);
			 return;
		 }

		 ndec = (int)tmpFloat1;
		 mc = (int)tmpFloat2;
		 count_read = count_read + 1;
		 t_f2 = post_write_uart(&token_succ,1);
		 //printf("Initialization successful!\n\r");
		 printf("ndec: %d, mc: %d, iniFlag: %f \n\r",ndec,mc,iniFlag);

		 if (iniFlag > 0.999 && iniFlag < 1.001) {
			 printf("H and A initialization. Memory reallocation.\n\r");

			 maxmem = ndec+1;
			 if (mc > ndec)
				 maxmem = mc+1;

			 if (hasCalloc == 1) {
				 // free the memory first
				 free(H); free(H_ori); free(invH); free(invH_ori); free(c); free(c_ori);
				 free(A); free(b); free(x_ini); free(Aw); free(bw); free(gx); free(lambda);
				 free(p); free(tmpMat1); free(tmpMat2); free(tmpMat3); free(tmpMat4); free(tmpMat5);
				 free(tmpMat6); free(tmpVec1); free(tmpVec2); free(iterPoint); free(xStar); free(x);
				 free(Af); free(Lv); free(Yv); free(Rv); free(pvStar); free(isInWl);
				 free(uv); free(vl); free(wv); free(lambdal); free(lambdaf); free(P); free(G);
				 free(wf); free(wl); free(orderPermu); free(lx); free(ux);
			 }

			 // Initialize QP parameters
			 H = (double*)calloc(ndec*ndec,sizeof(double));
			 H_ori = (double*)calloc(ndec*ndec,sizeof(double));
			 invH = (double*)calloc(ndec*ndec,sizeof(double));
			 invH_ori = (double*)calloc(ndec*ndec,sizeof(double));
			 c = (double*)calloc(ndec,sizeof(double));
			 c_ori = (double*)calloc(ndec,sizeof(double));
			 A = (double*)calloc(ndec*mc,sizeof(double));
			 b = (double*)calloc(mc,sizeof(double));
			 x_ini = (double*)calloc(ndec,sizeof(double));
			 Aw = (double*)calloc(ndec*(ndec+1),sizeof(double));
			 bw = (double*)calloc(ndec+1,sizeof(double));
			 gx = (double*)calloc(ndec+1,sizeof(double));
			 lambda = (double*)calloc(ndec+1,sizeof(double));
			 p = (double*)calloc(ndec+1,sizeof(double));
			 tmpMat1 = (double*)calloc((ndec+1)*ndec,sizeof(double));
			 tmpMat2 = (double*)calloc(ndec*ndec,sizeof(double));
			 tmpMat3 = (double*)calloc(ndec*ndec,sizeof(double));
			 tmpMat4 = (double*)calloc(ndec*ndec,sizeof(double));
			 tmpMat5 = (double*)calloc(ndec*ndec,sizeof(double));
			 tmpMat6 = (double*)calloc(ndec*ndec,sizeof(double));
			 tmpVec1 = (double*)calloc(maxmem,sizeof(double));
			 tmpVec2 = (double*)calloc(maxmem,sizeof(double));
			 iterPoint = (int*)calloc(1,sizeof(double));
			 xStar = (double*)calloc(ndec,sizeof(double));
			 x = (double*)calloc(ndec,sizeof(double));

			 Af = (double*)calloc(ndec*ndec,sizeof(double));
			 Lv = (double*)calloc(ndec*ndec,sizeof(double));
			 Yv = (double*)calloc(ndec*ndec,sizeof(double));
			 Rv = (double*)calloc(ndec*ndec,sizeof(double));
			 pvStar = (double*)calloc(ndec,sizeof(double));
			 isInWl = (int*)calloc(mc,sizeof(int));
			 uv = (double*)calloc(ndec,sizeof(double));
			 vl = (double*)calloc(ndec,sizeof(double));
			 wv = (double*)calloc(ndec,sizeof(double));
			 lambdal = (double*)calloc(ndec,sizeof(double));
			 lambdaf = (double*)calloc(ndec,sizeof(double));
			 P = (double*)calloc(ndec*ndec,sizeof(double));
			 G = (double*)calloc(4,sizeof(double));
			 wf = (int*)calloc(ndec,sizeof(int));
			 wl = (int*)calloc(ndec,sizeof(int));
			 orderPermu = (int*)calloc(ndec,sizeof(int));
			 lx = (double*)calloc((nbc+1)*1,sizeof(double));
			 ux = (double*)calloc((nbc+1)*1,sizeof(double));

			 hasCalloc = 1;

			 start_time = Timestamp_get32();
			 end_time = Timestamp_get32();
			 dead_time = end_time - start_time;

			 // Initialize H
			 //发送读数据指令
			 t_f2 = post_write_uart(&token_init,1);
			 for (i = 0; i < ndec*ndec; i++) {
				 //读8字节数据
				 t_f1 = post_read_uart(&Buffer_read[0],8,delay);
				 memcpy(&count_r_verify,&Buffer_read[0],sizeof(float));
				 memcpy(&tmpFloat1,&Buffer_read[4],sizeof(float));
				 if(count_r_verify!=count_read) {
					 printf("Count_read Verify failed! \n\r");
					 printf("Is read successful? %d \n\r",t_f1);
					 return;
				 }
				 H[i] = (double)tmpFloat1;
				 t_f2 = post_write_uart(&token_succ,1);
			 }
			 //printf("H is:\n\r"); show_matrix(H,ndec,ndec);printf("\n");
			 count_read = count_read + 1;

			 // Initialize c
			 //发送读数据指令
			 t_f2 = post_write_uart(&token_init,1);
			 for (i = 0; i < ndec; i++) {
				 //读8字节数据
				 t_f1 = post_read_uart(&Buffer_read[0],8,delay);
				 memcpy(&count_r_verify,&Buffer_read[0],sizeof(float));
				 memcpy(&tmpFloat1,&Buffer_read[4],sizeof(float));
				 if(count_r_verify!=count_read) {
					 printf("Count_read Verify failed! \n\r");
					 printf("Is read successful? %d \n\r",t_f1);
					 return;
				 }
				 c[i] = (double)tmpFloat1;
				 t_f2 = post_write_uart(&token_succ,1);
			 }
			 //printf("c is:\n\r"); show_matrix(c,ndec,1);printf("\n");
			 count_read = count_read + 1;

			 // Initialize A
			 //发送读数据指令
			 t_f2 = post_write_uart(&token_init,1);
			 for (i = 0; i < ndec*mc; i++) {
				 //读8字节数据
				 t_f1 = post_read_uart(&Buffer_read[0],8,delay);
				 memcpy(&count_r_verify,&Buffer_read[0],sizeof(float));
				 memcpy(&tmpFloat1,&Buffer_read[4],sizeof(float));
				 if(count_r_verify!=count_read) {
					 printf("Count_read Verify failed! \n\r");
					 printf("Is read successful? %d \n\r",t_f1);
					 return;
				 }
				 A[i] = (double)tmpFloat1;
				 t_f2 = post_write_uart(&token_succ,1);
			 }
			 //printf("A is:\n\r"); show_matrix(A,ndec,mc);printf("\n");
			 count_read = count_read + 1;

			 // Initialize b
			 //发送读数据指令
			 t_f2 = post_write_uart(&token_init,1);
			 for (i = 0; i < mc; i++) {
				 //读8字节数据
				 t_f1 = post_read_uart(&Buffer_read[0],8,delay);
				 memcpy(&count_r_verify,&Buffer_read[0],sizeof(float));
				 memcpy(&tmpFloat1,&Buffer_read[4],sizeof(float));
				 if(count_r_verify!=count_read) {
					 printf("Count_read Verify failed! \n\r");
					 printf("Is read successful? %d \n\r",t_f1);
					 return;
				 }
				 b[i] = (double)tmpFloat1;
				 t_f2 = post_write_uart(&token_succ,1);
			 }
			 //printf("b is:\n\r"); show_matrix(b,mc,1);printf("\n");
			 count_read = count_read + 1;

			 // Initialize invH
			 //发送读数据指令
			 t_f2 = post_write_uart(&token_init,1);
			 for (i = 0; i < ndec*ndec; i++) {
				 //读8字节数据
				 t_f1 = post_read_uart(&Buffer_read[0],8,delay);
				 memcpy(&count_r_verify,&Buffer_read[0],sizeof(float));
				 memcpy(&tmpFloat1,&Buffer_read[4],sizeof(float));
				 if(count_r_verify!=count_read) {
					 printf("Count_read Verify failed! \n\r");
					 printf("Is read successful? %d \n\r",t_f1);
					 return;
				 }
				 invH[i] = (double)tmpFloat1;
				 t_f2 = post_write_uart(&token_succ,1);
			 }
			 //printf("invH is:\n\r"); show_matrix(invH,ndec,ndec);printf("\n");
			 count_read = count_read + 1;
		 }
		 else if (iniFlag > 0.001 || iniFlag < -0.001) {
			 printf("Wrong iniFlag!\n\r");
			 return;
		 }
		 //==============================DSP读串口===============================
		 // Read c
		 t_f2 = post_write_uart(&token_read,1);
		 for (i = 0; i < ndec; i++) {
			 //读8字节数据
			 t_f1 = post_read_uart(&Buffer_read[0],8,delay);
			 memcpy(&count_r_verify,&Buffer_read[0],sizeof(float));
			 memcpy(&tmpFloat1,&Buffer_read[4],sizeof(float));
			 if(count_r_verify!=count_read) {
				 printf("count_r_verify is %f:! \n\r",count_r_verify);
				 printf("count_read is %f:! \n\r",count_read);
				 printf("Count_read Verify failed! \n\r");
				 printf("Is read successful? %d \n\r",t_f1);
				 return;
			 }
			 c[i] = (double)tmpFloat1;
			 t_f2 = post_write_uart(&token_succ,1);
		 }
		 //printf("c is:\n\r"); show_matrix(c,ndec,1);printf("\n");

		 // Read b
		 for (i = 0; i < mc; i++) {
			 //读8字节数据
			 t_f1 = post_read_uart(&Buffer_read[0],8,delay);
			 memcpy(&count_r_verify,&Buffer_read[0],sizeof(float));
			 memcpy(&tmpFloat1,&Buffer_read[4],sizeof(float));
			 if(count_r_verify!=count_read) {
				 printf("Count_read Verify failed! \n\r");
				 printf("Is read successful? %d \n\r",t_f1);
				 return;
			 }
			 b[i] = (double)tmpFloat1;
			 t_f2 = post_write_uart(&token_succ,1);
		 }
		 //printf("b is:\n\r"); show_matrix(b,mc,1);printf("\n");

		 // Read x
		 for (i = 0; i < ndec; i++) {
			 //读8字节数据
			 t_f1 = post_read_uart(&Buffer_read[0],8,delay);
			 memcpy(&count_r_verify,&Buffer_read[0],sizeof(float));
			 memcpy(&tmpFloat1,&Buffer_read[4],sizeof(float));
			 if(count_r_verify!=count_read) {
				 printf("Count_read Verify failed! \n\r");
				 printf("Is read successful? %d \n\r",t_f1);
				 return;
			 }
			 x[i] = (double)tmpFloat1;
			 t_f2 = post_write_uart(&token_succ,1);
		 }
		 //printf("x is:\n\r"); show_matrix(x,ndec,1);printf("\n");

		 // Read wSize and w
		 t_f1 = post_read_uart(&Buffer_read[0],8,delay);
		 memcpy(&count_r_verify,&Buffer_read[0],sizeof(float));
		 memcpy(&tmpFloat1,&Buffer_read[4],sizeof(float));
		 if(count_r_verify!=count_read) {
			 printf("Count_read Verify failed! \n\r");
			 printf("Is read successful? %d \n\r",t_f1);
			 return;
		 }
		 wSize = (int)tmpFloat1;
		 t_f2 = post_write_uart(&token_succ,1);

		 for (i = 0; i < wSize; i++) {
			 //读8字节数据
			 t_f1 = post_read_uart(&Buffer_read[0],8,delay);
			 memcpy(&count_r_verify,&Buffer_read[0],sizeof(float));
			 memcpy(&tmpFloat1,&Buffer_read[4],sizeof(float));
			 if(count_r_verify!=count_read) {
				 printf("Count_read Verify failed! \n\r");
				 printf("Is read successful? %d \n\r",t_f1);
				 return;
			 }
			 w[i] = (int)tmpFloat1;
			 t_f2 = post_write_uart(&token_succ,1);
		 }
		 //printf("w is:\n\r"); show_matrixInt(w,wSize,1);printf("\n");

		 // Read solverOpt
		 t_f1 = post_read_uart(&Buffer_read[0],8,delay);
		 memcpy(&count_r_verify,&Buffer_read[0],sizeof(float));
		 memcpy(&tmpFloat1,&Buffer_read[4],sizeof(float));
		 if(count_r_verify!=count_read) {
			 printf("Count_read Verify failed! \n\r");
			 printf("Is read successful? %d \n\r",t_f1);
			 return;
		 }
		 solverOpt = (int)tmpFloat1;
		 t_f2 = post_write_uart(&token_succ,1);
		 count_read = count_read + 1;

	    //==============================DSP读串口结束===============================

		 MPC_start_time = Timestamp_get32();

		 // 在这里调用 QP 求解器
		 if (solverOpt == 1) {
			 printf("Invoke ASM.\n\r");
			 ASM(H, invH, c, A, b,
					 x, ndec, mc, w, wSize, Aw, bw,
					 H_ori, invH_ori, c_ori, p, gx, lambda, tmpVec1, tmpVec2,
					 tmpMat1, tmpMat2, tmpMat3, tmpMat4,
					 xStar, iterPoint, maxIter);
		 }
		 else if (solverOpt == 2) {
			 // Note that here we always set an empty initial working set
			 printf("Invoke WGS.\n\r");
			 wgsQP(H, c, A, lx, ux, b,
					 wf, wl, 0, 0, x, ndec, 0, mc,
					 orderPermu, H_ori, c_ori, Af, Lv, Yv, Rv,
					 pvStar, p, gx, isInWl,uv,vl,wv,
					 lambdal,lambdaf,lambda, P,G,
					 tmpVec1, tmpVec2, tmpMat1, tmpMat2, tmpMat3, tmpMat4, tmpMat5, tmpMat6,
					 xStar,iterPoint,maxIter) ;
		 }
		 else {
			 printf("Wrong solver option!\n\r");
			 return;
		 }

		 MPC_end_time = Timestamp_get32();
		 MPC_total_time = MPC_end_time - MPC_start_time - dead_time;
		 Time_MPC = (double)(MPC_total_time);
		 printf("time is(ms): %f\n",Time_MPC/freq.lo*1000);
		 Time_MPC_Total += Time_MPC/freq.lo*1000;
		 if (Time_MPC/freq.lo*1000 > Time_MPC_Worst)
			 Time_MPC_Worst = Time_MPC/freq.lo*1000;

		 //printf("x is:\n\r"); show_matrix(x,ndec,1); printf("\n\r");
		 printf("iter is:%d\n\r",iterPoint[0]);

		 //Here the DSP communicate with Matlab to transmite u(k) and get y(k)

		 //==============================DSP写串口==================================
		 //发送写数据指令
		 t_f2 = post_write_uart(&token_write,1);
		 for (i = 0; i < ndec; i++) {
			 //封装
			 tmpFloat1 = (float)x[i];
			 memcpy(&Buffer_write[0],&count_write,sizeof(float));
			 memcpy(&Buffer_write[4],&tmpFloat1,sizeof(float));
			 //写12位数据
			 t_f2 = post_write_uart(&Buffer_write[0],8);
		 }
		 tmpFloat1 = (float)iterPoint[0];
		 memcpy(&Buffer_write[4],&tmpFloat1,sizeof(float));
		 t_f2 = post_write_uart(&Buffer_write[0],8);
		 tmpFloat1 = (float)(Time_MPC/freq.lo*1000);
		 memcpy(&Buffer_write[4],&tmpFloat1,sizeof(float));
		 t_f2 = post_write_uart(&Buffer_write[0],8);
		 // count_write = count_write + 1;
		 // One successful computation, restart
		 count_write = 1;
		 count_read = 1;


	 }


	 // Memory free

	 printf("exit taskFxn()\n");
 }

 Void main()
 {
 //    Task_Handle task;
     Error_Block eb;

     System_printf("enter main()\n");

     Error_init(&eb);
 //    task = Task_create(taskFxn, NULL, &eb);
 //    if (task == NULL) {
 //        System_printf("Task_create() failed!\n");
 //        BIOS_exit(0);
 //    }

     BIOS_start();     /* enable interrupts and start SYS/BIOS */
 }
