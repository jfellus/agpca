 /**
 * \file Library.java
 * \author Philippe H. Gosselin
 * \version 4.0
 */

package retin.toolbox.core;

import java.io.*;

public class Library 
{
	static	String	defaultMode = "release";
	static	String	ccName;
	static	String	archName;

	public static	void	setDefaultMode (String mode)
	{
		defaultMode = mode;
	}

	public static	void 	setDefaultMode (String[] args)
	{
		for (int i=0;i<args.length;i++)
		{
			int j = args[i].indexOf("--mode=");
			if (j == 0)
				defaultMode = args[i].substring(7);
		}
	}

	public static	void	load (String libname)
	{
		String error;
		error = load (libname,defaultMode);
		if (error == null)
			return;
		if (defaultMode != "release") {
			error = load (libname,"release");
			if (error == null)
				return;
		}
		if (defaultMode != "devel") {
			error = load (libname,"devel");
			if (error == null)
				return;
		}
		if (defaultMode != "debug") {
			error = load (libname,"debug");
			if (error == null)
				return;	
		}

		try {
			System.out.println ("Error loading library "+libname+" with arch="+getArchName()+" and cc="+getCCName());
			System.out.println (error);
		}
		catch (Exception ex) { 
			System.out.println(ex.getMessage());
		}
		System.exit(2);
	}

	public static	String	load (String libname,String mode)
	{
		try {
			String lib = "retin_"+libname+"-"+getArchName()+"-"+getCCName()+"-"+mode;
//			System.out.println ("Library "+lib+" try...");
			System.loadLibrary(lib);
			System.out.println ("Library "+lib+" loaded");
		}
		catch (Exception ex) { return ex.getMessage(); }
		catch (UnsatisfiedLinkError ex) { 
                        ex.printStackTrace();
                        return ex.getMessage();
                }
		return null;		
	}

	public static	String	getCCName () throws Exception
	{
		if (ccName == null)
		{
			Process process = Runtime.getRuntime().exec("gcc -dumpversion");
			process.waitFor ();
			if (process.exitValue() != 0)
				throw new Exception ("Error gcc -dumpversion");
		    	BufferedReader inStream = new BufferedReader(new InputStreamReader(process.getInputStream()));  
         		String[] list = inStream.readLine().split("\\.");
			ccName = "gcc"+list[0]+list[1]+"-mt";
			process.destroy ();
		}
		return ccName;
	}

	public static	String	getArchName () throws Exception
	{
		if (archName == null)
		{
			archName = "lin";

		// bits
			Process process = Runtime.getRuntime().exec("uname -m");
			process.waitFor ();
			if (process.exitValue() != 0)
				throw new Exception ("Error arch");
			BufferedReader inStream = new BufferedReader(new InputStreamReader(process.getInputStream()));  
			String arch = inStream.readLine();
			if (arch.equals("i686"))
				archName += "32";
			else
				archName += "64";
			process.destroy ();

		// SSE
			process = Runtime.getRuntime().exec("cat /proc/cpuinfo");
			process.waitFor ();
			if (process.exitValue() != 0)
				throw new Exception ("Error /proc/cpuinfo");
			inStream = new BufferedReader(new InputStreamReader(process.getInputStream()));  
			String line,cpuinfo = "";
			while ( (line = inStream.readLine()) != null )
				cpuinfo += line;

			String	bestSSE = "sse";
			String[] sses = { "sse2", "sse3", "ssse3", "sse4", "avx" };
			for (int i=0;i<sses.length;i++)
			{
				if (cpuinfo.indexOf(sses[i]) >= 0)
					bestSSE = sses[i];
			}
			archName += bestSSE;
			process.destroy ();
		}
		return archName;
	}
};

