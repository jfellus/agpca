package retin.toolbox.core;

import java.io.File;
import java.io.FileInputStream;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Properties;

public class Modules {

    public static void init() {
        Modules.init(new File("init.properties"));
    }

    public static void init(File initFile)
    {
        try {
            if (!initFile.exists())
                return;

            Properties init = new Properties();
            FileInputStream in = new FileInputStream(initFile);
            init.load(in);
            in.close();

            Iterator<Entry<Object,Object>> ite = init.entrySet().iterator();
            while (ite.hasNext()) {
                Entry<Object,Object> entry = ite.next();
                String key = (String)entry.getKey();
                //System.out.println("Init module "+key);
                try {
                    Class.forName(key);
                }
                catch (ClassNotFoundException ex) {
                    System.err.println("No class "+key);
                }
                catch (Exception ex) {
                    System.out.println("Error with "+key);
                    ex.printStackTrace(System.err);
                }
            }
        }
        catch(Exception ex) {
            ex.printStackTrace(System.err);
        }
    }

}
